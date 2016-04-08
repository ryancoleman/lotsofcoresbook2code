# -*- coding: utf-8 -*-
from __future__ import print_function, division
import functools
from math import pi

import numpy as np
import ase.units as units
from ase.utils.timing import timer

import gpaw.fftw as fftw
from gpaw.band_descriptor import BandDescriptor
from gpaw.blacs import BlacsGrid, BlacsDescriptor, Redistributor
from gpaw.density import Density
from gpaw.hs_operators import MatrixOperator
from gpaw.lfc import BaseLFC
from gpaw.lcao.overlap import fbt
from gpaw.hamiltonian import Hamiltonian
from gpaw.matrix_descriptor import MatrixDescriptor
from gpaw.spherical_harmonics import Y, nablarlYL
from gpaw.spline import Spline
from gpaw.utilities import unpack, _fact as fac
from gpaw.utilities.blas import rk, r2k, gemv, gemm, axpy
from gpaw.utilities.progressbar import ProgressBar
from gpaw.wavefunctions.fdpw import FDPWWaveFunctions


class PW:
    name = 'pw'
    
    def __init__(self, ecut=340, fftwflags=fftw.ESTIMATE, cell=None):
        """Plane-wave basis mode.

        ecut: float
            Plane-wave cutoff in eV.
        fftwflags: int
            Flags for making FFTW plan (default is ESTIMATE).
        cell: 3x3 ndarray
            Use this unit cell to chose the planewaves."""

        self.ecut = ecut / units.Hartree
        self.fftwflags = fftwflags

        if cell is None:
            self.cell_cv = None
        else:
            self.cell_cv = cell / units.Bohr

    def __call__(self, diagksl, orthoksl, initksl, gd, *args):
        if self.cell_cv is None:
            ecut = self.ecut
        else:
            volume = abs(np.linalg.det(gd.cell_cv))
            volume0 = abs(np.linalg.det(self.cell_cv))
            ecut = self.ecut * (volume0 / volume)**(2 / 3.0)

        wfs = PWWaveFunctions(ecut, self.fftwflags,
                              diagksl, orthoksl, initksl, gd, *args)
        return wfs

    def __eq__(self, other):
        return (isinstance(other, PW) and self.ecut == other.ecut)

    def __ne__(self, other):
        return not self == other


class PWDescriptor:
    ndim = 1  # all 3d G-vectors are stored in a 1d ndarray

    def __init__(self, ecut, gd, dtype=None, kd=None,
                 fftwflags=fftw.ESTIMATE):

        assert gd.pbc_c.all()
        assert gd.comm.size == 1

        self.ecut = ecut
        self.gd = gd
        self.fftwflags = fftwflags

        N_c = gd.N_c
        self.comm = gd.comm

        assert ((gd.h_cv**2).sum(1) <= 0.5 * pi**2 / ecut).all()

        if dtype is None:
            if kd is None or kd.gamma:
                dtype = float
            else:
                dtype = complex
        self.dtype = dtype

        if dtype == float:
            Nr_c = N_c.copy()
            Nr_c[2] = N_c[2] // 2 + 1
            i_Qc = np.indices(Nr_c).transpose((1, 2, 3, 0))
            i_Qc[..., :2] += N_c[:2] // 2
            i_Qc[..., :2] %= N_c[:2]
            i_Qc[..., :2] -= N_c[:2] // 2
            self.tmp_Q = fftw.empty(Nr_c, complex)
            self.tmp_R = self.tmp_Q.view(float)[:, :, :N_c[2]]
        else:
            i_Qc = np.indices(N_c).transpose((1, 2, 3, 0))
            i_Qc += N_c // 2
            i_Qc %= N_c
            i_Qc -= N_c // 2
            self.tmp_Q = fftw.empty(N_c, complex)
            self.tmp_R = self.tmp_Q

        self.nbytes = self.tmp_R.nbytes

        self.fftplan = fftw.FFTPlan(self.tmp_R, self.tmp_Q, -1, fftwflags)
        self.ifftplan = fftw.FFTPlan(self.tmp_Q, self.tmp_R, 1, fftwflags)

        # Calculate reciprocal lattice vectors:
        B_cv = 2.0 * pi * gd.icell_cv
        i_Qc.shape = (-1, 3)
        self.G_Qv = np.dot(i_Qc, B_cv)
        self.nbytes += self.G_Qv.nbytes

        self.kd = kd
        if kd is None:
            self.K_qv = np.zeros((1, 3))
        else:
            self.K_qv = np.dot(kd.ibzk_qc, B_cv)

        # Map from vectors inside sphere to fft grid:
        self.Q_qG = []
        self.G2_qG = []
        Q_Q = np.arange(len(i_Qc), dtype=np.int32)

        self.ngmin = 100000000
        self.ngmax = 0
        for q, K_v in enumerate(self.K_qv):
            G2_Q = ((self.G_Qv + K_v)**2).sum(axis=1)
            mask_Q = (G2_Q <= 2 * ecut)
            if self.dtype == float:
                mask_Q &= ((i_Qc[:, 2] > 0) |
                           (i_Qc[:, 1] > 0) |
                           ((i_Qc[:, 0] >= 0) & (i_Qc[:, 1] == 0)))
            Q_G = Q_Q[mask_Q]
            self.Q_qG.append(Q_G)
            self.G2_qG.append(G2_Q[Q_G])
            ng = len(Q_G)
            self.ngmin = min(ng, self.ngmin)
            self.ngmax = max(ng, self.ngmax)
            self.nbytes += Q_G.nbytes + self.G2_qG[q].nbytes

        if kd is not None:
            self.ngmin = kd.comm.min(self.ngmin)
            self.ngmax = kd.comm.max(self.ngmax)

        self.n_c = np.array([self.ngmax])  # used by hs_operators.py XXX

    def get_reciprocal_vectors(self, q=0, add_q=True):
        """ Returns reciprocal lattice vectors plus q, G + q,
        in xyz coordinates.
        """

        assert q < len(self.K_qv), ('Choose a q-index belonging to ' +
                                    'the irreducible Brillouin zone.')
        q_v = self.K_qv[q]

        if add_q:
            G_Gv = self.G_Qv[self.Q_qG[q]] + q_v
        else:
            G_Gv = self.G_Qv[self.Q_qG[q]]

        return G_Gv

    def __getstate__(self):
        return (self.ecut, self.gd, self.dtype, self.kd, self.fftwflags)

    def __setstate__(self, state):
        self.__init__(*state)

    def estimate_memory(self, mem):
        mem.subnode('Arrays', self.nbytes)

    def bytecount(self, dtype=float):
        return self.ngmax * 16

    def zeros(self, x=(), dtype=None, q=-1):
        a_xG = self.empty(x, dtype, q)
        a_xG.fill(0.0)
        return a_xG

    def empty(self, x=(), dtype=None, q=-1):
        if dtype is not None:
            assert dtype == self.dtype
        if isinstance(x, int):
            x = (x,)
        if q == -1:
            shape = x + (self.ngmax,)
        else:
            shape = x + self.Q_qG[q].shape
        return np.empty(shape, complex)

    def fft(self, f_R, q=-1, Q_G=None):
        """Fast Fourier transform.

        Returns c(G) for G<Gc::

                   __
                  \        -iG.R
            c(G) = ) f(R) e
                  /__
                   R
        """

        self.tmp_R[:] = f_R

        self.fftplan.execute()
        if Q_G is None:
            Q_G = self.Q_qG[q]
        return self.tmp_Q.ravel()[Q_G]

    def ifft(self, c_G, q=-1):
        """Inverse fast Fourier transform.

        Returns::

                      __
                   1 \        iG.R
            f(R) = -  ) c(G) e
                   N /__
                      G
        """

        self.tmp_Q[:] = 0.0
        self.tmp_Q.ravel()[self.Q_qG[q]] = c_G
        if self.dtype == float:
            t = self.tmp_Q[:, :, 0]
            n, m = self.gd.N_c[:2] // 2 - 1
            t[0, -m:] = t[0, m:0:-1].conj()
            t[n:0:-1, -m:] = t[-n:, m:0:-1].conj()
            t[-n:, -m:] = t[n:0:-1, m:0:-1].conj()
            t[-n:, 0] = t[n:0:-1, 0].conj()
        self.ifftplan.execute()
        return self.tmp_R * (1.0 / self.tmp_R.size)

    def integrate(self, a_xg, b_yg=None,
                  global_integral=True, hermitian=False,
                  _transposed_result=None):
        """Integrate function(s) over domain.

        a_xg: ndarray
            Function(s) to be integrated.
        b_yg: ndarray
            If present, integrate a_xg.conj() * b_yg.
        global_integral: bool
            If the array(s) are distributed over several domains, then the
            total sum will be returned.  To get the local contribution
            only, use global_integral=False.
        hermitian: bool
            Result is hermitian.
        _transposed_result: ndarray
            Long story.  Don't use this unless you are a method of the
            MatrixOperator class ..."""

        if b_yg is None:
            # Only one array:
            assert self.dtype == float
            return a_xg[..., 0].real * self.gd.dv

        A_xg = a_xg.reshape((-1, a_xg.shape[-1]))
        B_yg = b_yg.reshape((-1, b_yg.shape[-1]))

        alpha = self.gd.dv / self.gd.N_c.prod()

        if self.dtype == float:
            alpha *= 2
            A_xg = A_xg.view(float)
            B_yg = B_yg.view(float)

        if _transposed_result is None:
            result_yx = np.zeros((len(B_yg), len(A_xg)), self.dtype)
        else:
            result_yx = _transposed_result

        if a_xg is b_yg:
            rk(alpha, A_xg, 0.0, result_yx)
        elif hermitian:
            r2k(0.5 * alpha, A_xg, B_yg, 0.0, result_yx)
        else:
            gemm(alpha, A_xg, B_yg, 0.0, result_yx, 'c')

        if self.dtype == float:
            correction_yx = np.outer(B_yg[:, 0], A_xg[:, 0])
            if hermitian:
                result_yx -= 0.25 * alpha * (correction_yx + correction_yx.T)
            else:
                result_yx -= 0.5 * alpha * correction_yx

        xshape = a_xg.shape[:-1]
        yshape = b_yg.shape[:-1]
        result = result_yx.T.reshape(xshape + yshape)

        if result.ndim == 0:
            return result.item()
        else:
            return result

    def interpolate(self, a_R, pd, q=-1):
        a_Q = self.tmp_Q
        b_Q = pd.tmp_Q

        e0, e1, e2 = 1 - self.gd.N_c % 2  # even or odd size
        a0, a1, a2 = pd.gd.N_c // 2 - self.gd.N_c // 2
        b0, b1, b2 = self.gd.N_c + (a0, a1, a2)

        if self.dtype == float:
            b2 = (b2 - a2) // 2 + 1
            a2 = 0
            axes = (0, 1)
        else:
            axes = (0, 1, 2)

        self.tmp_R[:] = a_R
        self.fftplan.execute()
        b_Q[:] = 0.0
        b_Q[a0:b0, a1:b1, a2:b2] = np.fft.fftshift(a_Q, axes=axes)

        if e0:
            b_Q[a0, a1:b1, a2:b2] *= 0.5
            b_Q[b0, a1:b1, a2:b2] = b_Q[a0, a1:b1, a2:b2]
            b0 += 1
        if e1:
            b_Q[a0:b0, a1, a2:b2] *= 0.5
            b_Q[a0:b0, b1, a2:b2] = b_Q[a0:b0, a1, a2:b2]
            b1 += 1
        if self.dtype == complex:
            if e2:
                b_Q[a0:b0, a1:b1, a2] *= 0.5
                b_Q[a0:b0, a1:b1, b2] = b_Q[a0:b0, a1:b1, a2]
        else:
            if e2:
                b_Q[a0:b0, a1:b1, b2 - 1] *= 0.5

        b_Q[:] = np.fft.ifftshift(b_Q, axes=axes)
        pd.ifftplan.execute()
        return pd.tmp_R * (1.0 / self.tmp_R.size), a_Q.ravel()[self.Q_qG[q]]

    def restrict(self, a_R, pd, q=-1):
        a_Q = pd.tmp_Q
        b_Q = self.tmp_Q

        e0, e1, e2 = 1 - pd.gd.N_c % 2  # even or odd size
        a0, a1, a2 = self.gd.N_c // 2 - pd.gd.N_c // 2
        b0, b1, b2 = pd.gd.N_c // 2 + self.gd.N_c // 2 + 1

        if self.dtype == float:
            b2 = pd.gd.N_c[2] // 2 + 1
            a2 = 0
            axes = (0, 1)
        else:
            axes = (0, 1, 2)

        self.tmp_R[:] = a_R
        self.fftplan.execute()
        b_Q[:] = np.fft.fftshift(b_Q, axes=axes)

        if e0:
            b_Q[a0, a1:b1, a2:b2] += b_Q[b0 - 1, a1:b1, a2:b2]
            b_Q[a0, a1:b1, a2:b2] *= 0.5
            b0 -= 1
        if e1:
            b_Q[a0:b0, a1, a2:b2] += b_Q[a0:b0, b1 - 1, a2:b2]
            b_Q[a0:b0, a1, a2:b2] *= 0.5
            b1 -= 1
        if self.dtype == complex and e2:
            b_Q[a0:b0, a1:b1, a2] += b_Q[a0:b0, a1:b1, b2 - 1]
            b_Q[a0:b0, a1:b1, a2] *= 0.5
            b2 -= 1

        a_Q[:] = b_Q[a0:b0, a1:b1, a2:b2]
        a_Q[:] = np.fft.ifftshift(a_Q, axes=axes)
        a_G = a_Q.ravel()[pd.Q_qG[q]] / 8
        pd.ifftplan.execute()
        return pd.tmp_R * (1.0 / self.tmp_R.size), a_G

    def map(self, pd, q=-1):
        N_c = np.array(self.tmp_Q.shape)
        N3_c = pd.tmp_Q.shape
        Q2_G = self.Q_qG[q]
        Q2_Gc = np.empty((len(Q2_G), 3), int)
        Q2_Gc[:, 0], r_G = divmod(Q2_G, N_c[1] * N_c[2])
        Q2_Gc.T[1:] = divmod(r_G, N_c[2])
        if self.dtype == float:
            C = 2
        else:
            C = 3
        Q2_Gc[:, :C] += N_c[:C] // 2
        Q2_Gc[:, :C] %= N_c[:C]
        Q2_Gc[:, :C] -= N_c[:C] // 2
        Q2_Gc[:, :C] %= N3_c[:C]
        Q3_G = Q2_Gc[:, 2] + N3_c[2] * (Q2_Gc[:, 1] + N3_c[1] * Q2_Gc[:, 0])
        G3_Q = np.empty(N3_c, int).ravel()
        G3_Q[pd.Q_qG[q]] = np.arange(len(pd.Q_qG[q]))
        return G3_Q[Q3_G]

    def gemm(self, alpha, psit_nG, C_mn, beta, newpsit_mG):
        """Helper function for MatrixOperator class."""
        if self.dtype == float:
            psit_nG = psit_nG.view(float)
            newpsit_mG = newpsit_mG.view(float)
        gemm(alpha, psit_nG, C_mn, beta, newpsit_mG)

    def gemv(self, alpha, psit_nG, C_n, beta, newpsit_G, trans='t'):
        """Helper function for CG eigensolver."""
        if self.dtype == float:
            psit_nG = psit_nG.view(float)
            newpsit_G = newpsit_G.view(float)
        gemv(alpha, psit_nG, C_n, beta, newpsit_G, trans)


def count_reciprocal_vectors(ecut, gd, q_c):
    N_c = gd.N_c
    i_Qc = np.indices(N_c).transpose((1, 2, 3, 0))
    i_Qc += N_c // 2
    i_Qc %= N_c
    i_Qc -= N_c // 2

    B_cv = 2.0 * pi * gd.icell_cv
    i_Qc.shape = (-1, 3)
    Gpq_Qv = np.dot(i_Qc, B_cv) + np.dot(q_c, B_cv)
    
    G2_Q = (Gpq_Qv**2).sum(axis=1)
    return (G2_Q <= 2 * ecut).sum()

            
class Preconditioner:
    """Preconditioner for KS equation.

    From:

      Teter, Payne and Allen, Phys. Rev. B 40, 12255 (1989)

    as modified by:

      Kresse and Furthmüller, Phys. Rev. B 54, 11169 (1996)
    """

    def __init__(self, G2_qG, pd):
        self.G2_qG = G2_qG
        self.pd = pd

    def calculate_kinetic_energy(self, psit_xG, kpt):
        G2_G = self.G2_qG[kpt.q]
        return [self.pd.integrate(0.5 * G2_G * psit_G, psit_G)
                for psit_G in psit_xG]

    def __call__(self, R_xG, kpt, ekin_x):
        G2_G = self.G2_qG[kpt.q]
        PR_xG = np.empty_like(R_xG)
        for PR_G, R_G, ekin in zip(PR_xG, R_xG, ekin_x):
            x_G = 1 / ekin / 3 * G2_G
            a_G = 27.0 + x_G * (18.0 + x_G * (12.0 + x_G * 8.0))
            PR_G[:] = 4.0 / 3 / ekin * R_G * a_G / (a_G + 16.0 * x_G**4)
        return -PR_xG


class PWWaveFunctions(FDPWWaveFunctions):
    mode = 'pw'

    def __init__(self, ecut, fftwflags,
                 diagksl, orthoksl, initksl,
                 gd, nvalence, setups, bd, dtype,
                 world, kd, kptband_comm, timer):
        self.ecut = ecut
        self.fftwflags = fftwflags

        self.ng_k = None  # number of G-vectors for all IBZ k-points

        FDPWWaveFunctions.__init__(self, diagksl, orthoksl, initksl,
                                   gd, nvalence, setups, bd, dtype,
                                   world, kd, kptband_comm, timer)

        self.orthoksl.gd = self.pd
        self.matrixoperator = MatrixOperator(self.orthoksl)

    def empty(self, n=(), global_array=False, realspace=False,
              q=-1):
        if realspace:
            return self.gd.empty(n, self.dtype, global_array)
        else:
            return self.pd.empty(n, self.dtype, q)

    def integrate(self, a_xg, b_yg=None, global_integral=True):
        return self.pd.integrate(a_xg, b_yg, global_integral)

    def bytes_per_wave_function(self):
        return 16 * self.pd.ngmax

    def set_setups(self, setups):
        self.timer.start('PWDescriptor')
        self.pd = PWDescriptor(self.ecut, self.gd, self.dtype, self.kd,
                               self.fftwflags)
        self.timer.stop('PWDescriptor')

        # Build array of number of plane wave coefficiants for all k-points
        # in the IBZ:
        self.ng_k = np.zeros(self.kd.nibzkpts, dtype=int)
        for kpt in self.kpt_u:
            if kpt.s == 0:
                self.ng_k[kpt.k] = len(self.pd.Q_qG[kpt.q])
        self.kd.comm.sum(self.ng_k)

        self.pt = PWLFC([setup.pt_j for setup in setups], self.pd)

        FDPWWaveFunctions.set_setups(self, setups)

    def summary(self, fd):
        fd.write('Wave functions: Plane wave expansion\n')
        fd.write('      Cutoff energy: %.3f eV\n' %
                 (self.pd.ecut * units.Hartree))
        if self.dtype == float:
            fd.write('      Number of coefficients: %d (reduced to %d)\n' %
                     (self.pd.ngmax * 2 - 1, self.pd.ngmax))
        else:
            fd.write('      Number of coefficients (min, max): %d, %d\n' %
                     (self.pd.ngmin, self.pd.ngmax))
        if fftw.FFTPlan is fftw.NumpyFFTPlan:
            fd.write("      Using Numpy's FFT\n")
        else:
            fd.write('      Using FFTW library\n')

    def make_preconditioner(self, block=1):
        return Preconditioner(self.pd.G2_qG, self.pd)

    def apply_pseudo_hamiltonian(self, kpt, hamiltonian, psit_xG, Htpsit_xG):
        """Apply the non-pseudo Hamiltonian i.e. without PAW corrections."""
        Htpsit_xG[:] = 0.5 * self.pd.G2_qG[kpt.q] * psit_xG
        for psit_G, Htpsit_G in zip(psit_xG, Htpsit_xG):
            psit_R = self.pd.ifft(psit_G, kpt.q)
            Htpsit_G += self.pd.fft(psit_R * hamiltonian.vt_sG[kpt.s], kpt.q)

    def add_orbital_density(self, nt_G, kpt, n):
        axpy(1.0, abs(self.pd.ifft(kpt.psit_nG[n], kpt.q))**2, nt_G)

    def add_to_density_from_k_point_with_occupation(self, nt_sR, kpt, f_n):
        nt_R = nt_sR[kpt.s]
        for f, psit_G in zip(f_n, kpt.psit_nG):
            nt_R += f * abs(self.pd.ifft(psit_G, kpt.q))**2

    def calculate_kinetic_energy_density(self):
        if self.kpt_u[0].f_n is None:
            raise RuntimeError

        taut_sR = self.gd.zeros(self.nspins)
        for kpt in self.kpt_u:
            G_Gv = self.pd.get_reciprocal_vectors(q=kpt.q)
            for f, psit_G in zip(kpt.f_n, kpt.psit_nG):
                for v in range(3):
                    taut_sR[kpt.s] += 0.5 * f * abs(
                        self.pd.ifft(1j * G_Gv[:, v] * psit_G, kpt.q))**2

        self.kd.comm.sum(taut_sR)
        self.band_comm.sum(taut_sR)
        return taut_sR

    def apply_mgga_orbital_dependent_hamiltonian(self, kpt, psit_xG,
                                                 Htpsit_xG, dH_asp,
                                                 dedtaut_R):
        G_Gv = self.pd.get_reciprocal_vectors(q=kpt.q)
        for psit_G, Htpsit_G in zip(psit_xG, Htpsit_xG):
            for v in range(3):
                a_R = self.pd.ifft(1j * G_Gv[:, v] * psit_G, kpt.q)
                axpy(-0.5, 1j * G_Gv[:, v] *
                     self.pd.fft(dedtaut_R * a_R, kpt.q),
                     Htpsit_G)

    def _get_wave_function_array(self, u, n, realspace=True, phase=None):
        psit_G = FDPWWaveFunctions._get_wave_function_array(self, u, n,
                                                            realspace)
        if not realspace:
            zeropadded_G = np.zeros(self.pd.ngmax, complex)
            zeropadded_G[:len(psit_G)] = psit_G
            return zeropadded_G

        kpt = self.kpt_u[u]
        if self.kd.gamma:
            return self.pd.ifft(psit_G)
        else:
            if phase is None:
                N_c = self.gd.N_c
                k_c = self.kd.ibzk_kc[kpt.k]
                eikr_R = np.exp(2j * pi * np.dot(np.indices(N_c).T,
                                                 k_c / N_c).T)
            else:
                eikr_R = phase
            return self.pd.ifft(psit_G, kpt.q) * eikr_R

    def get_wave_function_array(self, n, k, s, realspace=True,
                                cut=True):
        psit_G = FDPWWaveFunctions.get_wave_function_array(self, n, k, s,
                                                           realspace)
        if cut and psit_G is not None and not realspace:
            psit_G = psit_G[:self.ng_k[k]].copy()

        return psit_G

    def write(self, writer, write_wave_functions=False):
        writer['Mode'] = 'pw'
        writer['PlaneWaveCutoff'] = self.ecut

        if not write_wave_functions:
            return

        writer.dimension('nplanewaves', self.pd.ngmax)
        writer.add('PseudoWaveFunctions',
                   ('nspins', 'nibzkpts', 'nbands', 'nplanewaves'),
                   dtype=complex)

        for s in range(self.nspins):
            for k in range(self.kd.nibzkpts):
                for n in range(self.bd.nbands):
                    psit_G = self.get_wave_function_array(n, k, s,
                                                          realspace=False,
                                                          cut=False)
                    writer.fill(psit_G, s, k, n)

        writer.add('PlaneWaveIndices', ('nibzkpts', 'nplanewaves'),
                   dtype=np.int32)

        if self.bd.comm.rank > 0:
            return

        Q_G = np.empty(self.pd.ngmax, np.int32)
        for r in range(self.kd.comm.size):
            for q, ks in enumerate(self.kd.get_indices(r)):
                s, k = divmod(ks, self.kd.nibzkpts)
                ng = self.ng_k[k]
                if s == 1:
                    return
                if r == self.kd.comm.rank:
                    Q_G[:ng] = self.pd.Q_qG[q]
                    if r > 0:
                        self.kd.comm.send(Q_G, 0)
                if self.kd.comm.rank == 0:
                    if r > 0:
                        self.kd.comm.receive(Q_G, r)
                    Q_G[ng:] = -1
                    writer.fill(Q_G, k)

    def read(self, reader, hdf5):
        assert reader['version'] >= 3
        for kpt in self.kpt_u:
            if kpt.s == 0:
                Q_G = reader.get('PlaneWaveIndices', kpt.k)
                ng = self.ng_k[kpt.k]
                assert (Q_G[:ng] == self.pd.Q_qG[kpt.q]).all()
                assert (Q_G[ng:] == -1).all()

        assert not hdf5
        if self.bd.comm.size == 1:
            for kpt in self.kpt_u:
                ng = self.ng_k[kpt.k]
                kpt.psit_nG = reader.get_reference('PseudoWaveFunctions',
                                                   (kpt.s, kpt.k),
                                                   length=ng)
            return

        for kpt in self.kpt_u:
            kpt.psit_nG = self.empty(self.bd.mynbands, q=kpt.q)
            ng = self.ng_k[kpt.k]
            for myn, psit_G in enumerate(kpt.psit_nG):
                n = self.bd.global_index(myn)
                psit_G[:] = reader.get('PseudoWaveFunctions',
                                       kpt.s, kpt.k, n)[..., :ng]

    def hs(self, ham, q=-1, s=0, md=None):
        npw = len(self.pd.Q_qG[q])
        N = self.pd.tmp_R.size

        if md is None:
            H_GG = np.zeros((npw, npw), complex)
            S_GG = np.zeros((npw, npw), complex)
            G1 = 0
            G2 = npw
        else:
            H_GG = md.zeros(dtype=complex)
            S_GG = md.zeros(dtype=complex)
            G1, G2 = md.my_blocks(S_GG).next()[:2]

        H_GG.ravel()[G1::npw + 1] = (0.5 * self.pd.gd.dv / N *
                                     self.pd.G2_qG[q][G1:G2])

        for G in range(G1, G2):
            x_G = self.pd.zeros(q=q)
            x_G[G] = 1.0
            H_GG[G - G1] += (self.pd.gd.dv / N *
                             self.pd.fft(ham.vt_sG[s] *
                                         self.pd.ifft(x_G, q), q))

        S_GG.ravel()[G1::npw + 1] = self.pd.gd.dv / N

        f_IG = self.pt.expand(q)
        nI = len(f_IG)
        dH_II = np.zeros((nI, nI))
        dS_II = np.zeros((nI, nI))
        I1 = 0
        for a in self.pt.my_atom_indices:
            dH_ii = unpack(ham.dH_asp[a][s])
            dS_ii = self.setups[a].dO_ii
            I2 = I1 + len(dS_ii)
            dH_II[I1:I2, I1:I2] = dH_ii / N**2
            dS_II[I1:I2, I1:I2] = dS_ii / N**2
            I1 = I2

        H_GG += np.dot(f_IG.T[G1:G2].conj(), np.dot(dH_II, f_IG))
        S_GG += np.dot(f_IG.T[G1:G2].conj(), np.dot(dS_II, f_IG))

        return H_GG, S_GG

    @timer('Full diag')
    def diagonalize_full_hamiltonian(self, ham, atoms, occupations, txt,
                                     nbands=None, scalapack=None, expert=False):
        assert self.dtype == complex

        if nbands is None:
            nbands = self.pd.ngmin // self.bd.comm.size * self.bd.comm.size
        else:
            assert nbands <= self.pd.ngmin

        if expert:
            iu = nbands
        else:
            iu = None

        self.bd = bd = BandDescriptor(nbands, self.bd.comm)

        p = functools.partial(print, file=txt)
        p('Diagonalizing full Hamiltonian ({0} lowest bands)'.format(nbands))
        p('Matrix size (min, max): {0}, {1}'.format(self.pd.ngmin,
                                                    self.pd.ngmax))
        mem = 3 * self.pd.ngmax**2 * 16 / bd.comm.size / 1024**2
        p('Approximate memory usage per core: {0:.3f} MB'.format(mem))
        if bd.comm.size > 1:
            if isinstance(scalapack, (list, tuple)):
                nprow, npcol, b = scalapack
            else:
                nprow = int(round(bd.comm.size**0.5))
                while bd.comm.size % nprow != 0:
                    nprow -= 1
                npcol = bd.comm.size // nprow
                b = 64
            p('ScaLapack grid: {0}x{1},'.format(nprow, npcol),
              'block-size:', b)
            bg = BlacsGrid(bd.comm, bd.comm.size, 1)
            bg2 = BlacsGrid(bd.comm, nprow, npcol)
            scalapack = True
        else:
            nprow = npcol = 1
            scalapack = False

        self.pt.set_positions(atoms.get_scaled_positions())
        self.kpt_u[0].P_ani = None
        self.allocate_arrays_for_projections(self.pt.my_atom_indices)

        myslice = bd.get_slice()

        pb = ProgressBar(txt)
        nkpt = len(self.kpt_u)
        
        for u, kpt in enumerate(self.kpt_u):
            pb.update(u / nkpt)
            npw = len(self.pd.Q_qG[kpt.q])
            if scalapack:
                mynpw = -(-npw // bd.comm.size)
                md = BlacsDescriptor(bg, npw, npw, mynpw, npw)
                md2 = BlacsDescriptor(bg2, npw, npw, b, b)
            else:
                md = md2 = MatrixDescriptor(npw, npw)

            with self.timer('Build H and S'):
                H_GG, S_GG = self.hs(ham, kpt.q, kpt.s, md)

            if scalapack:
                r = Redistributor(bd.comm, md, md2)
                H_GG = r.redistribute(H_GG)
                S_GG = r.redistribute(S_GG)

            psit_nG = md2.empty(dtype=complex)
            eps_n = np.empty(npw)

            with self.timer('Diagonalize'):
                if not scalapack:
                    md2.general_diagonalize_dc(H_GG, S_GG, psit_nG, eps_n, 
                                               iu=iu)
                else:
                    md2.general_diagonalize_dc(H_GG, S_GG, psit_nG, eps_n)
            del H_GG, S_GG

            kpt.eps_n = eps_n[myslice].copy()

            if scalapack:
                md3 = BlacsDescriptor(bg, npw, npw, bd.mynbands, npw)
                r = Redistributor(bd.comm, md2, md3)
                psit_nG = r.redistribute(psit_nG)

            kpt.psit_nG = psit_nG[:bd.mynbands].copy()
            del psit_nG

            with self.timer('Projections'):
                self.pt.integrate(kpt.psit_nG, kpt.P_ani, kpt.q)

            kpt.f_n = None

        pb.finish()
        
        occupations.calculate(self)

    def initialize_from_lcao_coefficients(self, basis_functions, mynbands):
        N_c = self.gd.N_c

        psit_nR = self.gd.empty(mynbands, self.dtype)

        for kpt in self.kpt_u:
            if self.kd.gamma:
                emikr_R = 1.0
            else:
                k_c = self.kd.ibzk_kc[kpt.k]
                emikr_R = np.exp(-2j * pi *
                                 np.dot(np.indices(N_c).T, k_c / N_c).T)

            psit_nR[:] = 0.0
            basis_functions.lcao_to_grid(kpt.C_nM, psit_nR, kpt.q)
            kpt.C_nM = None

            kpt.psit_nG = self.pd.empty(self.bd.mynbands, q=kpt.q)
            for n in range(mynbands):
                kpt.psit_nG[n] = self.pd.fft(psit_nR[n] * emikr_R, kpt.q)

    def random_wave_functions(self, mynao):
        rs = np.random.RandomState(self.world.rank)
        for kpt in self.kpt_u:
            psit_nG = kpt.psit_nG[mynao:]
            weight_G = 1.0 / (1.0 + self.pd.G2_qG[kpt.q])
            psit_nG.real = rs.uniform(-1, 1, psit_nG.shape) * weight_G
            psit_nG.imag = rs.uniform(-1, 1, psit_nG.shape) * weight_G

    def estimate_memory(self, mem):
        FDPWWaveFunctions.estimate_memory(self, mem)
        self.pd.estimate_memory(mem.subnode('PW-descriptor'))

    def get_kinetic_stress(self):
        sigma_vv = np.zeros((3, 3), dtype=complex)
        pd = self.pd
        dOmega = pd.gd.dv / pd.gd.N_c.prod()
        if pd.dtype == float:
            dOmega *= 2
        K_qv = self.pd.K_qv
        for kpt in self.kpt_u:
            G_Gv = pd.get_reciprocal_vectors(q=kpt.q, add_q=False)
            psit2_G = 0.0
            for n, f in enumerate(kpt.f_n):
                psit2_G += f * np.abs(kpt.psit_nG[n])**2
            for alpha in range(3):
                Ga_G = G_Gv[:, alpha] + K_qv[kpt.q, alpha]
                for beta in range(3):
                    Gb_G = G_Gv[:, beta] + K_qv[kpt.q, beta]
                    sigma_vv[alpha, beta] += (psit2_G * Ga_G * Gb_G).sum()

        sigma_vv *= -dOmega

        self.bd.comm.sum(sigma_vv)
        self.kd.comm.sum(sigma_vv)
        return sigma_vv


def ft(spline):
    l = spline.get_angular_momentum_number()
    rc = 50.0
    N = 2**10
    assert spline.get_cutoff() <= rc

    dr = rc / N
    r_r = np.arange(N) * dr
    dk = pi / 2 / rc
    k_q = np.arange(2 * N) * dk
    f_r = spline.map(r_r) * (4 * pi)

    f_q = fbt(l, f_r, r_r, k_q)
    f_q[1:] /= k_q[1:]**(2 * l + 1)
    f_q[0] = (np.dot(f_r, r_r**(2 + 2 * l)) *
              dr * 2**l * fac[l] / fac[2 * l + 1])

    return Spline(l, k_q[-1], f_q)


class PWLFC(BaseLFC):
    def __init__(self, spline_aj, pd, blocksize=5000):
        """Reciprocal-space plane-wave localized function collection.

        spline_aj: list of list of spline objects
            Splines.
        pd: PWDescriptor
            Plane-wave descriptor object.
        blocksize: int
            Block-size to use when looping over G-vectors.  Default is to do
            all G-vectors in one big block."""

        self.pd = pd
        self.spline_aj = spline_aj

        self.dtype = pd.dtype

        self.initialized = False
        
        # These will be filled in later:
        self.lf_aj = []
        self.Y_qLG = []

        if blocksize is not None:
            if pd.ngmax <= blocksize:
                # No need to block G-vectors
                blocksize = None
        self.blocksize = blocksize

        # These are set later in set_potitions():
        self.eikR_qa = None
        self.my_atom_indices = None
        self.indices = None
        self.pos_av = None
        self.nI = None

    def initialize(self):
        if self.initialized:
            return
            
        cache = {}
        lmax = -1

        # Fourier transform radial functions:
        for a, spline_j in enumerate(self.spline_aj):
            self.lf_aj.append([])
            for spline in spline_j:
                l = spline.get_angular_momentum_number()
                if spline not in cache:
                    f = ft(spline)
                    f_qG = []
                    for G2_G in self.pd.G2_qG:
                        G_G = G2_G**0.5
                        f_qG.append(f.map(G_G))
                    cache[spline] = f_qG
                else:
                    f_qG = cache[spline]
                self.lf_aj[a].append((l, f_qG))
                lmax = max(lmax, l)

        # Spherical harmonics:
        for q, K_v in enumerate(self.pd.K_qv):
            G_Gv = self.pd.get_reciprocal_vectors(q=q)
            Y_LG = np.empty(((lmax + 1)**2, len(G_Gv)))
            for L in range((lmax + 1)**2):
                Y_LG[L] = Y(L, *G_Gv.T)
            self.Y_qLG.append(Y_LG)
            
        self.initialized = True

    def estimate_memory(self, mem):
        splines = set()
        lmax = -1
        for spline_j in self.spline_aj:
            for spline in spline_j:
                splines.add(spline)
                l = spline.get_angular_momentum_number()
                lmax = max(lmax, l)
        nbytes = ((len(splines) + (lmax + 1)**2) *
                  sum(G2_G.nbytes for G2_G in self.pd.G2_qG))
        mem.subnode('Arrays', nbytes)

    def get_function_count(self, a):
        return sum(2 * l + 1 for l, f_qG in self.lf_aj[a])

    def __iter__(self):
        I1 = 0
        for a in self.my_atom_indices:
            j = 0
            i1 = 0
            for l, f_qG in self.lf_aj[a]:
                i2 = i1 + 2 * l + 1
                I2 = I1 + 2 * l + 1
                yield a, j, i1, i2, I1, I2
                i1 = i2
                I1 = I2
                j += 1

    def set_positions(self, spos_ac):
        self.initialize()
        kd = self.pd.kd
        if kd is None or kd.gamma:
            self.eikR_qa = np.ones((1, len(spos_ac)))
        else:
            self.eikR_qa = np.exp(2j * pi * np.dot(kd.ibzk_qc, spos_ac.T))

        self.pos_av = np.dot(spos_ac, self.pd.gd.cell_cv)

        self.my_atom_indices = np.arange(len(spos_ac))
        self.indices = []
        I1 = 0
        for a in self.my_atom_indices:
            I2 = I1 + self.get_function_count(a)
            self.indices.append((a, I1, I2))
            I1 = I2
        self.nI = I1

    def expand(self, q=-1, G1=0, G2=None):
        if G2 is None:
            G2 = self.Y_qLG[q].shape[1]
        f_IG = np.empty((self.nI, G2 - G1), complex)
        emiGR_Ga = np.exp(-1j * np.dot(self.pd.G_Qv[self.pd.Q_qG[q][G1:G2]],
                                       self.pos_av.T))
        for a, j, i1, i2, I1, I2 in self:
            l, f_qG = self.lf_aj[a][j]
            f_IG[I1:I2] = (emiGR_Ga[:, a] * f_qG[q][G1:G2] * (-1.0j)**l *
                           self.Y_qLG[q][l**2:(l + 1)**2, G1:G2])
        return f_IG

    def block(self, q=-1):
        nG = self.Y_qLG[q].shape[1]
        if self.blocksize:
            G1 = 0
            while G1 < nG:
                G2 = min(G1 + self.blocksize, nG)
                yield G1, G2
                G1 = G2
        else:
            yield 0, nG

    def add(self, a_xG, c_axi=1.0, q=-1, f0_IG=None):
        if isinstance(c_axi, float):
            assert q == -1, a_xG.dims == 1
            a_xG += (c_axi / self.pd.gd.dv) * self.expand(-1).sum(0)
            return

        c_xI = np.empty(a_xG.shape[:-1] + (self.nI,), self.pd.dtype)
        for a, I1, I2 in self.indices:
            c_xI[..., I1:I2] = c_axi[a] * self.eikR_qa[q][a].conj()
        c_xI = c_xI.reshape((np.prod(c_xI.shape[:-1], dtype=int), self.nI))

        a_xG = a_xG.reshape((-1, a_xG.shape[-1])).view(self.pd.dtype)

        for G1, G2 in self.block(q):
            if f0_IG is None:
                f_IG = self.expand(q, G1, G2)
            else:
                f_IG = f0_IG

            if self.pd.dtype == float:
                f_IG = f_IG.view(float)
                G1 *= 2
                G2 *= 2

            gemm(1.0 / self.pd.gd.dv, f_IG, c_xI, 1.0, a_xG[:, G1:G2])

    def integrate(self, a_xG, c_axi=None, q=-1):
        c_xI = np.zeros(a_xG.shape[:-1] + (self.nI,), self.pd.dtype)

        b_xI = c_xI.reshape((np.prod(c_xI.shape[:-1], dtype=int), self.nI))
        a_xG = a_xG.reshape((-1, a_xG.shape[-1]))

        alpha = 1.0 / self.pd.gd.N_c.prod()
        if self.pd.dtype == float:
            alpha *= 2
            a_xG = a_xG.view(float)

        if c_axi is None:
            c_axi = self.dict(a_xG.shape[:-1])

        x = 0.0
        for G1, G2 in self.block(q):
            f_IG = self.expand(q, G1, G2)
            if self.pd.dtype == float:
                if G1 == 0:
                    f_IG[:, 0] *= 0.5
                f_IG = f_IG.view(float)
                G1 *= 2
                G2 *= 2

            gemm(alpha, f_IG, a_xG[:, G1:G2], x, b_xI, 'c')
            x = 1.0

        for a, I1, I2 in self.indices:
            c_axi[a][:] = self.eikR_qa[q][a] * c_xI[..., I1:I2]

        return c_axi

    def derivative(self, a_xG, c_axiv, q=-1):
        c_vxI = np.zeros((3,) + a_xG.shape[:-1] + (self.nI,), self.pd.dtype)
        b_vxI = c_vxI.reshape((3, np.prod(c_vxI.shape[1:-1], dtype=int),
                               self.nI))
        a_xG = a_xG.reshape((-1, a_xG.shape[-1])).view(self.pd.dtype)

        alpha = 1.0 / self.pd.gd.N_c.prod()

        K_v = self.pd.K_qv[q]

        x = 0.0
        for G1, G2 in self.block(q):
            f_IG = self.expand(q, G1, G2)
            G_Gv = self.pd.G_Qv[self.pd.Q_qG[q][G1:G2]]
            if self.pd.dtype == float:
                for v in range(3):
                    gemm(2 * alpha,
                         (f_IG * 1.0j * G_Gv[:, v]).view(float),
                         a_xG[:, 2 * G1:2 * G2],
                         x, b_vxI[v], 'c')
            else:
                for v in range(3):
                    gemm(-alpha,
                         f_IG * (G_Gv[:, v] + K_v[v]),
                         a_xG[:, G1:G2],
                         x, b_vxI[v], 'c')
            x = 1.0

        for v in range(3):
            if self.pd.dtype == float:
                for a, I1, I2 in self.indices:
                    c_axiv[a][..., v] = c_vxI[v, ..., I1:I2]
            else:
                for a, I1, I2 in self.indices:
                    c_axiv[a][..., v] = (1.0j * self.eikR_qa[q][a] *
                                         c_vxI[v, ..., I1:I2])

    def stress_tensor_contribution(self, a_xG, c_axi=1.0, q=-1):
        cache = {}
        for a, j, i1, i2, I1, I2 in self:
            spline = self.spline_aj[a][j]
            if spline not in cache:
                s = ft(spline)
                G_G = self.pd.G2_qG[q]**0.5
                f_G = []
                dfdGoG_G = []
                for G in G_G:
                    f, dfdG = s.get_value_and_derivative(G)
                    if G < 1e-10:
                        G = 1.0
                    f_G.append(f)
                    dfdGoG_G.append(dfdG / G)
                f_G = np.array(f_G)
                dfdGoG_G = np.array(dfdGoG_G)
                cache[spline] = (f_G, dfdGoG_G)
            else:
                f_G, dfdGoG_G = cache[spline]

        if isinstance(c_axi, float):
            c_axi = dict((a, c_axi) for a in range(len(self.pos_av)))

        G0_Gv = self.pd.get_reciprocal_vectors(q=q)

        stress_vv = np.zeros((3, 3))
        for G1, G2 in self.block(q):
            G_Gv = G0_Gv[G1:G2]
            aa_xG = a_xG[..., G1:G2]
            for v1 in range(3):
                for v2 in range(3):
                    stress_vv[v1, v2] += self._stress_tensor_contribution(
                        v1, v2, cache, G1, G2, G_Gv, aa_xG, c_axi, q)

        return stress_vv

    def _stress_tensor_contribution(self, v1, v2, cache, G1, G2,
                                    G_Gv, a_xG, c_axi, q):
        f_IG = np.empty((self.nI, G2 - G1), complex)
        K_v = self.pd.K_qv[q]
        for a, j, i1, i2, I1, I2 in self:
            l = self.lf_aj[a][j][0]
            spline = self.spline_aj[a][j]
            f_G, dfdGoG_G = cache[spline]

            emiGR_G = np.exp(-1j * np.dot(G_Gv, self.pos_av[a]))
            f_IG[I1:I2] = (emiGR_G * (-1.0j)**l *
                           np.exp(1j * np.dot(K_v, self.pos_av[a])) * (
                               dfdGoG_G[G1:G2] * G_Gv[:, v1] * G_Gv[:, v2] *
                               self.Y_qLG[q][l**2:(l + 1)**2, G1:G2] +
                               f_G[G1:G2] * G_Gv[:, v1] *
                               [nablarlYL(L, G_Gv.T)[v2]
                                for L in range(l**2, (l + 1)**2)]))

        c_xI = np.zeros(a_xG.shape[:-1] + (self.nI,), self.pd.dtype)

        b_xI = c_xI.reshape((np.prod(c_xI.shape[:-1], dtype=int), self.nI))
        a_xG = a_xG.reshape((-1, a_xG.shape[-1]))

        alpha = 1.0 / self.pd.gd.N_c.prod()
        if self.pd.dtype == float:
            alpha *= 2
            if G1 == 0:
                f_IG[:, 0] *= 0.5
            f_IG = f_IG.view(float)
            a_xG = a_xG.view(float)

        gemm(alpha, f_IG, a_xG, 0.0, b_xI, 'c')

        stress = 0.0
        for a, I1, I2 in self.indices:
            stress -= self.eikR_qa[q][a] * (c_axi[a] * c_xI[..., I1:I2]).sum()
        return stress.real


class PseudoCoreKineticEnergyDensityLFC(PWLFC):
    def add(self, tauct_R):
        tauct_R += self.pd.ifft(1.0 / self.pd.gd.dv * self.expand(-1).sum(0))

    def derivative(self, dedtaut_R, dF_aiv):
        PWLFC.derivative(self, self.pd.fft(dedtaut_R), dF_aiv)


class ReciprocalSpaceDensity(Density):
    def __init__(self, gd, finegd, nspins, charge, collinear=True):
        Density.__init__(self, gd, finegd, nspins, charge, collinear)

        self.ecut2 = 0.5 * pi**2 / (self.gd.h_cv**2).sum(1).max() * 0.9999
        self.pd2 = PWDescriptor(self.ecut2, self.gd)
        self.ecut3 = 0.5 * pi**2 / (self.finegd.h_cv**2).sum(1).max() * 0.9999
        self.pd3 = PWDescriptor(self.ecut3, self.finegd)

        self.G3_G = self.pd2.map(self.pd3)

    def initialize(self, setups, timer, magmom_av, hund):
        Density.initialize(self, setups, timer, magmom_av, hund)

        spline_aj = []
        for setup in setups:
            if setup.nct is None:
                spline_aj.append([])
            else:
                spline_aj.append([setup.nct])
        self.nct = PWLFC(spline_aj, self.pd2)

        self.ghat = PWLFC([setup.ghat_l for setup in setups], self.pd3)

    def set_positions(self, spos_ac, rank_a=None):
        Density.set_positions(self, spos_ac, rank_a)
        self.nct_q = self.pd2.zeros()
        self.nct.add(self.nct_q, 1.0 / self.nspins)
        self.nct_G = self.pd2.ifft(self.nct_q)

    def interpolate_pseudo_density(self, comp_charge=None):
        """Interpolate pseudo density to fine grid."""
        if comp_charge is None:
            comp_charge = self.calculate_multipole_moments()

        if self.nt_sg is None:
            self.nt_sg = self.finegd.empty(self.nspins * self.ncomp**2)
            self.nt_sQ = self.pd2.empty(self.nspins * self.ncomp**2)

        for nt_G, nt_Q, nt_g in zip(self.nt_sG, self.nt_sQ, self.nt_sg):
            nt_g[:], nt_Q[:] = self.pd2.interpolate(nt_G, self.pd3)

    def interpolate(self, in_xR, out_xR=None):
        """Interpolate array(s)."""
        if out_xR is None:
            out_xR = self.finegd.empty(in_xR.shape[:-3])

        a_xR = in_xR.reshape((-1,) + in_xR.shape[-3:])
        b_xR = out_xR.reshape((-1,) + out_xR.shape[-3:])

        for in_R, out_R in zip(a_xR, b_xR):
            out_R[:] = self.pd2.interpolate(in_R, self.pd3)[0]

        return out_xR

    def calculate_pseudo_charge(self):
        self.nt_Q = self.nt_sQ[:self.nspins].sum(axis=0)
        self.rhot_q = self.pd3.zeros()
        self.rhot_q[self.G3_G] = self.nt_Q * 8
        self.ghat.add(self.rhot_q, self.Q_aL)
        self.rhot_q[0] = 0.0

    def get_pseudo_core_kinetic_energy_density_lfc(self):
        return PseudoCoreKineticEnergyDensityLFC(
            [[setup.tauct] for setup in self.setups], self.pd2)

    def calculate_dipole_moment(self):
        if np.__version__ < '1.6.0':
            raise NotImplementedError
        pd = self.pd3
        N_c = pd.tmp_Q.shape

        m0_q, m1_q, m2_q = [i_G == 0
                            for i_G in np.unravel_index(pd.Q_qG[0], N_c)]
        rhot_q = self.rhot_q.imag
        rhot_cs = [rhot_q[m1_q & m2_q],
                   rhot_q[m0_q & m2_q],
                   rhot_q[m0_q & m1_q]]
        d_c = [np.dot(rhot_s[1:], 1.0 / np.arange(1, len(rhot_s)))
               for rhot_s in rhot_cs]
        return -np.dot(d_c, pd.gd.cell_cv) / pi * pd.gd.dv


class ReciprocalSpaceHamiltonian(Hamiltonian):
    def __init__(self, gd, finegd, pd2, pd3, nspins, setups, timer, xc,
                 world, kptband_comm, vext=None, collinear=True):
        Hamiltonian.__init__(self, gd, finegd, nspins, setups, timer, xc,
                             world, kptband_comm, vext, collinear)

        self.vbar = PWLFC([[setup.vbar] for setup in setups], pd2)
        self.pd2 = pd2
        self.pd3 = pd3

        class PS:
            def initialize(self):
                pass

            def get_stencil(self):
                return '????'

            def estimate_memory(self, mem):
                pass

        self.poisson = PS()
        self.npoisson = 0

    def summary(self, fd):
        Hamiltonian.summary(self, fd)
        fd.write('Interpolation: FFT\n')
        fd.write('Poisson solver: FFT\n')

    def set_positions(self, spos_ac, rank_a=None):
        Hamiltonian.set_positions(self, spos_ac, rank_a)
        self.vbar_Q = self.pd2.zeros()
        self.vbar.add(self.vbar_Q)

    def update_pseudo_potential(self, density):
        self.ebar = self.pd2.integrate(self.vbar_Q, density.nt_sQ.sum(0))

        self.timer.start('Poisson')
        self.vHt_q = 4 * pi * density.rhot_q
        self.vHt_q[1:] /= self.pd3.G2_qG[0][1:]
        self.epot = 0.5 * self.pd3.integrate(self.vHt_q, density.rhot_q)
        self.timer.stop('Poisson')

        # Calculate atomic hamiltonians:
        W_aL = {}
        for a in density.D_asp:
            W_aL[a] = np.empty((self.setups[a].lmax + 1)**2)
        density.ghat.integrate(self.vHt_q, W_aL)

        self.vt_Q = self.vbar_Q + self.vHt_q[density.G3_G] / 8
        self.vt_sG[:] = self.pd2.ifft(self.vt_Q)

        self.timer.start('XC 3D grid')
        vxct_sg = self.finegd.zeros(self.nspins)
        self.exc = self.xc.calculate(self.finegd, density.nt_sg, vxct_sg)

        for vt_G, vxct_g in zip(self.vt_sG, vxct_sg):
            vxc_G, vxc_Q = self.pd3.restrict(vxct_g, self.pd2)
            vt_G += vxc_G
            self.vt_Q += vxc_Q / self.nspins
        self.timer.stop('XC 3D grid')

        ekin = 0.0
        for vt_G, nt_G in zip(self.vt_sG, density.nt_sG):
            ekin -= self.gd.integrate(vt_G, nt_G)
        ekin += self.gd.integrate(self.vt_sG, density.nct_G).sum()

        eext = 0.0

        return ekin, self.epot, self.ebar, eext, self.exc, W_aL

    def restrict(self, in_xR, out_xR=None):
        """Restrict array."""
        if out_xR is None:
            out_xR = self.gd.empty(in_xR.shape[:-3])

        a_xR = in_xR.reshape((-1,) + in_xR.shape[-3:])
        b_xR = out_xR.reshape((-1,) + out_xR.shape[-3:])

        for in_R, out_R in zip(a_xR, b_xR):
            out_R[:] = self.pd3.restrict(in_R, self.pd2)[0]

        return out_xR

    def calculate_forces2(self, dens, ghat_aLv, nct_av, vbar_av):
        dens.ghat.derivative(self.vHt_q, ghat_aLv)
        dens.nct.derivative(self.vt_Q, nct_av)
        self.vbar.derivative(dens.nt_sQ.sum(0), vbar_av)
