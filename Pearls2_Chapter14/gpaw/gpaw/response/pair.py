from __future__ import print_function
import sys
from math import pi

import numpy as np
from ase.units import Hartree
from ase.utils import devnull
from ase.utils.timing import timer, Timer

import gpaw.mpi as mpi
from gpaw import GPAW
from gpaw.fd_operators import Gradient
from gpaw.occupations import FermiDirac
from gpaw.response.math_func import two_phi_planewave_integrals
from gpaw.utilities.blas import gemm
from gpaw.wavefunctions.pw import PWLFC
import gpaw.io.tar as io

import warnings


class KPoint:
    def __init__(self, s, K, n1, n2, blocksize, na, nb,
                 ut_nR, eps_n, f_n, P_ani, shift_c):
        self.s = s    # spin index
        self.K = K    # BZ k-point index
        self.n1 = n1  # first band
        self.n2 = n2  # first band not included
        self.blocksize = blocksize
        self.na = na  # first band of block
        self.nb = nb  # first band of block not included
        self.ut_nR = ut_nR      # periodic part of wave functions in real-space
        self.eps_n = eps_n      # eigenvalues
        self.f_n = f_n          # occupation numbers
        self.P_ani = P_ani      # PAW projections
        self.shift_c = shift_c  # long story - see the
        # PairDensity.construct_symmetry_operators() method


class PairDensity:
    def __init__(self, calc, ecut=50,
                 ftol=1e-6, threshold=1,
                 real_space_derivatives=False,
                 world=mpi.world, txt=sys.stdout, timer=None, nblocks=1,
                 gate_voltage=None):
        if ecut is not None:
            ecut /= Hartree

        if gate_voltage is not None:
            gate_voltage /= Hartree

        self.ecut = ecut
        self.ftol = ftol
        self.threshold = threshold
        self.real_space_derivatives = real_space_derivatives
        self.world = world
        self.gate_voltage = gate_voltage

        if nblocks == 1:
            self.blockcomm = self.world.new_communicator([world.rank])
            self.kncomm = world
        else:
            assert world.size % nblocks == 0, world.size
            rank1 = world.rank // nblocks * nblocks
            rank2 = rank1 + nblocks
            self.blockcomm = self.world.new_communicator(range(rank1, rank2))
            ranks = range(world.rank % nblocks, world.size, nblocks)
            self.kncomm = self.world.new_communicator(ranks)

        if world.rank != 0:
            txt = devnull
        elif isinstance(txt, str):
            txt = open(txt, 'w')
        self.fd = txt

        self.timer = timer or Timer()

        if isinstance(calc, str):
            print('Reading ground state calculation:\n  %s' % calc,
                  file=self.fd)
            if not calc.split('.')[-1] == 'gpw':
                calc = calc + '.gpw'
            self.reader = io.Reader(calc, comm=mpi.serial_comm)
            calc = GPAW(calc, txt=None, communicator=mpi.serial_comm,
                        read_projections=False)
        else:
            self.reader = None
            assert calc.wfs.world.size == 1

        assert calc.wfs.kd.symmetry.symmorphic
        self.calc = calc

        if gate_voltage is not None:
            self.add_gate_voltage(gate_voltage)

        self.spos_ac = calc.atoms.get_scaled_positions()

        self.nocc1 = None  # number of completely filled bands
        self.nocc2 = None  # number of non-empty bands
        self.count_occupied_bands()

        self.vol = abs(np.linalg.det(calc.wfs.gd.cell_cv))

        self.ut_sKnvR = None  # gradient of wave functions for optical limit

        print('Number of blocks:', nblocks, file=self.fd)

    def add_gate_voltage(self, gate_voltage=0):
        """Shifts the Fermi-level by e * Vg. By definition e = 1."""
        assert isinstance(self.calc.occupations, FermiDirac)
        print('Shifting Fermi-level by %.2f eV' % (gate_voltage * Hartree),
              file=self.fd)

        for kpt in self.calc.wfs.kpt_u:
            kpt.f_n = (self.shift_occupations(kpt.eps_n, gate_voltage)
                       * kpt.weight)

    def shift_occupations(self, eps_n, gate_voltage):
        """Shift fermilevel."""
        fermi = self.calc.occupations.get_fermi_level() + gate_voltage
        width = self.calc.occupations.width
        tmp = (eps_n - fermi) / width
        f_n = np.zeros_like(eps_n)
        f_n[tmp <= 100] = 1 / (1 + np.exp(tmp[tmp <= 100]))
        f_n[tmp > 100] = 0.0
        return f_n

    def count_occupied_bands(self):
        self.nocc1 = 9999999
        self.nocc2 = 0
        for kpt in self.calc.wfs.kpt_u:
            f_n = kpt.f_n / kpt.weight
            self.nocc1 = min((f_n > 1 - self.ftol).sum(), self.nocc1)
            self.nocc2 = max((f_n > self.ftol).sum(), self.nocc2)
        print('Number of completely filled bands:', self.nocc1, file=self.fd)
        print('Number of partially filled bands:', self.nocc2, file=self.fd)
        print('Total number of bands:', self.calc.wfs.bd.nbands,
              file=self.fd)

    def distribute_k_points_and_bands(self, band1, band2, kpts=None):
        """Distribute spins, k-points and bands.

        nbands: int
            Number of bands for each spin/k-point combination.

        The attribute self.mysKn1n2 will be set to a list of (s, K, n1, n2)
        tuples that this process handles.
        """

        wfs = self.calc.wfs

        if kpts is None:
            kpts = range(wfs.kd.nbzkpts)

        nbands = band2 - band1
        size = self.kncomm.size
        rank = self.kncomm.rank
        ns = wfs.nspins
        nk = len(kpts)
        n = (ns * nk * nbands + size - 1) // size
        i1 = rank * n
        i2 = min(i1 + n, ns * nk * nbands)

        self.mysKn1n2 = []
        i = 0
        for s in range(ns):
            for K in kpts:
                n1 = min(max(0, i1 - i), nbands)
                n2 = min(max(0, i2 - i), nbands)
                if n1 != n2:
                    self.mysKn1n2.append((s, K, n1 + band1, n2 + band1))
                i += nbands

        print('BZ k-points:', self.calc.wfs.kd.description, file=self.fd)
        print('Distributing spins, k-points and bands (%d x %d x %d)' %
              (ns, nk, nbands),
              'over %d process%s' %
              (self.kncomm.size, ['es', ''][self.kncomm.size == 1]),
              file=self.fd)
        print('Number of blocks:', self.blockcomm.size, file=self.fd)

    @timer('Get a k-point')
    def get_k_point(self, s, K, n1, n2, block=False):
        """Return wave functions for a specific k-point and spin.

        s: int
            Spin index (0 or 1).
        K: int
            BZ k-point index.
        n1, n2: int
            Range of bands to include.
        """
        wfs = self.calc.wfs

        if block:
            nblocks = self.blockcomm.size
            rank = self.blockcomm.rank
        else:
            nblocks = 1
            rank = 0
        blocksize = (n2 - n1 + nblocks - 1) // nblocks
        na = min(n1 + rank * blocksize, n2)
        nb = min(na + blocksize, n2)

        U_cc, T, a_a, U_aii, shift_c, time_reversal = \
            self.construct_symmetry_operators(K)
        ik = wfs.kd.bz2ibz_k[K]
        kpt = wfs.kpt_u[s * wfs.kd.nibzkpts + ik]

        eps_n = kpt.eps_n[n1:n2]
        f_n = kpt.f_n[n1:n2] / kpt.weight

        psit_nG = kpt.psit_nG
        ut_nR = wfs.gd.empty(nb - na, wfs.dtype)
        for n in range(na, nb):
            ut_nR[n - na] = T(wfs.pd.ifft(psit_nG[n], ik))
            
        P_ani = []
        i1 = 0
        if self.reader is not None:
            P_nI = self.reader.get('Projections', kpt.s, kpt.k)
        else:
            iP_ani = kpt.P_ani

        for b, U_ii in zip(a_a, U_aii):
            i2 = i1 + len(U_ii)
            if self.reader is not None:
                P_ni = np.dot(P_nI[na:nb, i1:i2], U_ii)
            else:
                P_ni = np.dot(iP_ani[b][na:nb], U_ii)
            if time_reversal:
                P_ni = P_ni.conj()
            P_ani.append(P_ni)
            i1 = i2
        
        return KPoint(s, K, n1, n2, blocksize, na, nb,
                      ut_nR, eps_n, f_n, P_ani, shift_c)

    @timer('Calculate pair-densities')
    def calculate_pair_densities(self, ut1cc_R, C1_aGi, kpt2, pd, Q_G):
        """Calculate FFT of pair-densities and add PAW corrections.

        ut1cc_R: 3-d complex ndarray
            Complex conjugate of the periodic part of the left hand side
            wave function.
        C1_aGi: list of ndarrays
            PAW corrections for all atoms.
        kpt2: KPoint object
            Right hand side k-point object.
        pd: PWDescriptor
            Plane-wave descriptor for for q=k2-k1.
        Q_G: 1-d int ndarray
            Mapping from flattened 3-d FFT grid to 0.5(G+q)^2<ecut sphere.
        """

        dv = pd.gd.dv
        n_mG = pd.empty(kpt2.blocksize)
        myblocksize = kpt2.nb - kpt2.na

        for ut_R, n_G in zip(kpt2.ut_nR, n_mG):
            n_R = ut1cc_R * ut_R
            with self.timer('fft'):
                n_G[:] = pd.fft(n_R, 0, Q_G) * dv

        # PAW corrections:
        with self.timer('gemm'):
            for C1_Gi, P2_mi in zip(C1_aGi, kpt2.P_ani):
                gemm(1.0, C1_Gi, P2_mi, 1.0, n_mG[:myblocksize], 't')

        if self.blockcomm.size == 1:
            return n_mG
        else:
            n_MG = pd.empty(kpt2.blocksize * self.blockcomm.size)
            self.blockcomm.all_gather(n_mG, n_MG)
            return n_MG[:kpt2.n2 - kpt2.n1]

    def get_fft_indices(self, K1, K2, q_c, pd, shift0_c):
        """Get indices for G-vectors inside cutoff sphere."""
        kd = self.calc.wfs.kd
        N_G = pd.Q_qG[0]
        shift_c = (shift0_c +
                   (q_c - kd.bzk_kc[K2] + kd.bzk_kc[K1]).round().astype(int))
        if shift_c.any():
            n_cG = np.unravel_index(N_G, pd.gd.N_c)
            n_cG = [n_G + shift for n_G, shift in zip(n_cG, shift_c)]
            N_G = np.ravel_multi_index(n_cG, pd.gd.N_c, 'wrap')
        return N_G

    def construct_symmetry_operators(self, K):
        """Construct symmetry operators for wave function and PAW projections.

        We want to transform a k-point in the irreducible part of the BZ to
        the corresponding k-point with index K.

        Returns U_cc, T, a_a, U_aii, shift_c and time_reversal, where:

        * U_cc is a rotation matrix.
        * T() is a function that transforms the periodic part of the wave
          function.
        * a_a is a list of symmetry related atom indices
        * U_aii is a list of rotation matrices for the PAW projections
        * shift_c is three integers: see code below.
        * time_reversal is a flag - if True, projections should be complex
          conjugated.

        See the get_k_point() method for how to use these tuples.
        """

        wfs = self.calc.wfs
        kd = wfs.kd

        s = kd.sym_k[K]
        U_cc = kd.symmetry.op_scc[s]
        time_reversal = kd.time_reversal_k[K]
        ik = kd.bz2ibz_k[K]
        k_c = kd.bzk_kc[K]
        ik_c = kd.ibzk_kc[ik]

        sign = 1 - 2 * time_reversal
        shift_c = np.dot(U_cc, ik_c) - k_c * sign
        assert np.allclose(shift_c.round(), shift_c)
        shift_c = shift_c.round().astype(int)

        if (U_cc == np.eye(3)).all():
            T = lambda f_R: f_R
        else:
            N_c = self.calc.wfs.gd.N_c
            i_cr = np.dot(U_cc.T, np.indices(N_c).reshape((3, -1)))
            i = np.ravel_multi_index(i_cr, N_c, 'wrap')
            T = lambda f_R: f_R.ravel()[i].reshape(N_c)

        if time_reversal:
            T0 = T
            T = lambda f_R: T0(f_R).conj()
            shift_c *= -1

        a_a = []
        U_aii = []
        for a, id in enumerate(self.calc.wfs.setups.id_a):
            b = kd.symmetry.a_sa[s, a]
            S_c = np.dot(self.spos_ac[a], U_cc) - self.spos_ac[b]
            x = np.exp(2j * pi * np.dot(ik_c, S_c))
            U_ii = wfs.setups[a].R_sii[s].T * x
            a_a.append(b)
            U_aii.append(U_ii)

        return U_cc, T, a_a, U_aii, shift_c, time_reversal

    @timer('Initialize PAW corrections')
    def initialize_paw_corrections(self, pd, soft=False):
        wfs = self.calc.wfs
        q_v = pd.K_qv[0]
        optical_limit = np.allclose(q_v, 0)

        G_Gv = pd.get_reciprocal_vectors()
        if optical_limit:
            G_Gv[0] = 1

        pos_av = np.dot(self.spos_ac, pd.gd.cell_cv)

        # Collect integrals for all species:
        Q_xGii = {}
        for id, atomdata in wfs.setups.setups.items():
            if soft:
                ghat = PWLFC([atomdata.ghat_l], pd)
                ghat.set_positions(np.zeros((1, 3)))
                Q_LG = ghat.expand()
                Q_Gii = np.dot(atomdata.Delta_iiL, Q_LG).T
            else:
                Q_Gii = two_phi_planewave_integrals(G_Gv, atomdata)
                ni = atomdata.ni
                Q_Gii.shape = (-1, ni, ni)

            Q_xGii[id] = Q_Gii

        Q_aGii = []
        for a, atomdata in enumerate(wfs.setups):
            id = wfs.setups.id_a[a]
            Q_Gii = Q_xGii[id]
            x_G = np.exp(-1j * np.dot(G_Gv, pos_av[a]))
            Q_aGii.append(x_G[:, np.newaxis, np.newaxis] * Q_Gii)
            if optical_limit:
                Q_aGii[a][0] = atomdata.dO_ii

        return Q_aGii

    @timer('Optical limit')
    def update_optical_limit(self, n, m, kpt1, kpt2, deps_m, df_m, n_mG):
        if self.ut_sKnvR is None or kpt1.K not in self.ut_sKnvR[kpt1.s]:
            self.ut_sKnvR = self.calculate_derivatives(kpt1)

        # Relative threshold for perturbation theory
        threshold = self.threshold

        kd = self.calc.wfs.kd
        gd = self.calc.wfs.gd
        k_c = kd.bzk_kc[kpt1.K] + kpt1.shift_c
        k_v = 2 * np.pi * np.dot(k_c, np.linalg.inv(gd.cell_cv).T)

        ut_vR = self.ut_sKnvR[kpt1.s][kpt1.K][n]
        atomdata_a = self.calc.wfs.setups
        C_avi = [np.dot(atomdata.nabla_iiv.T, P_ni[n])
                 for atomdata, P_ni in zip(atomdata_a, kpt1.P_ani)]

        blockbands = kpt2.nb - kpt2.na
        n0_mv = np.empty((kpt2.blocksize, 3), dtype=complex)
        nt_m = np.empty(kpt2.blocksize, dtype=complex)
        n0_mv[:blockbands] = -self.calc.wfs.gd.integrate(ut_vR, kpt2.ut_nR).T
        nt_m[:blockbands] = self.calc.wfs.gd.integrate(kpt1.ut_nR[n],
                                                       kpt2.ut_nR)
        n0_mv += 1j * nt_m[:, np.newaxis] * k_v[np.newaxis, :]

        for C_vi, P_mi in zip(C_avi, kpt2.P_ani):
            P_mi = P_mi.copy()
            gemm(1.0, C_vi, P_mi, 1.0, n0_mv[:blockbands], 'c')

        if self.blockcomm.size != 1:
            n0_Mv = np.empty((kpt2.blocksize * self.blockcomm.size, 3),
                             dtype=complex)
            self.blockcomm.all_gather(n0_mv, n0_Mv)
            n0_mv = n0_Mv[:kpt2.n2 - kpt2.n1]

        # In case not all unoccupied bands are included
        if n0_mv.shape[0] != len(m):
            n0_mv = n0_mv[m]

        deps_m = deps_m.copy()
        deps_m[deps_m >= 0.0] = np.inf

        smallness_mv = np.abs(-1e-3 * n0_mv / deps_m[:, np.newaxis])
        inds_mv = (np.logical_and(np.inf > smallness_mv,
                                  smallness_mv > threshold))

        if inds_mv.any():
            indent8 = ' ' * 8
            print('\n    WARNING: Optical limit perturbation' +
                  ' theory failed for:', file=self.fd)
            print(indent8 + 'kpt_c = [%1.2f, %1.2f, %1.2f]'
                  % (k_c[0], k_c[1], k_c[2]), file=self.fd)
            inds_m = inds_mv.any(axis=1)
            depsi_m = deps_m[inds_m]
            n0i_mv = np.abs(n0_mv[inds_m])
            smallness_mv = smallness_mv[inds_m]
            for depsi, n0i_v, smallness_v in zip(depsi_m, n0i_mv,
                                                 smallness_mv):
                print(indent8 + 'Energy eigenvalue difference %1.2e ' % -depsi,
                      file=self.fd)
                print(indent8 + 'Matrix element' +
                      ' %1.2e %1.2e %1.2e' % (n0i_v[0], n0i_v[1], n0i_v[2]),
                      file=self.fd)
                print(indent8 + 'Smallness' +
                      ' %1.2e %1.2e %1.2e\n' % (smallness_v[0],
                                                smallness_v[1],
                                                smallness_v[2]),
                      file=self.fd)

        n0_mv *= 1j / deps_m[:, np.newaxis]
        n0_mv[inds_mv] = 0
        n_mG[:, 0] = n0_mv[:, 0]

        return n0_mv

    @timer('Intraband')
    def update_intraband(self, ind_m, kpt):
        kd = self.calc.wfs.kd
        gd = self.calc.wfs.gd
        k_c = kd.bzk_kc[kpt.K] + kpt.shift_c
        k_v = 2 * np.pi * np.dot(k_c, np.linalg.inv(gd.cell_cv).T)
        atomdata_a = self.calc.wfs.setups

        ut_mvR = self.calc.wfs.gd.zeros((len(ind_m), 3), complex)
        for ind, ut_vR in zip(ind_m, ut_mvR):
            ut_vR[:] = self.make_derivative(kpt.s, kpt.K,
                                            kpt.na + ind,
                                            kpt.na + ind + 1)[0]
        npartocc = len(ind_m)
        ut_mR = kpt.ut_nR[ind_m]

        nabla0_mmv = np.zeros((npartocc, npartocc, 3), dtype=complex)
        for m in range(npartocc):
            ut_vR = ut_mvR[m]
            C_avi = [np.dot(atomdata.nabla_iiv.T, P_mi[ind_m[m]])
                     for atomdata, P_mi in zip(atomdata_a, kpt.P_ani)]

            nabla0_mv = -self.calc.wfs.gd.integrate(ut_vR, ut_mR).T
            nt_m = self.calc.wfs.gd.integrate(ut_mR[m], ut_mR)
            nabla0_mv += 1j * nt_m[:, np.newaxis] * k_v[np.newaxis, :]

            for C_vi, P_mi in zip(C_avi, kpt.P_ani):
                gemm(1.0, C_vi, P_mi[ind_m[0:npartocc]], 1.0, nabla0_mv, 'c')

            nabla0_mmv[m] = nabla0_mv

        return nabla0_mmv

    def calculate_derivatives(self, kpt):
        ut_sKnvR = [{}, {}]
        ut_nvR = self.make_derivative(kpt.s, kpt.K, kpt.n1, kpt.n2)
        ut_sKnvR[kpt.s][kpt.K] = ut_nvR

        return ut_sKnvR

    @timer('Derivatives')
    def make_derivative(self, s, K, n1, n2):
        wfs = self.calc.wfs
        if self.real_space_derivatives:
            grad_v = [Gradient(wfs.gd, v, 1.0, 4, complex).apply
                      for v in range(3)]

        U_cc, T, a_a, U_aii, shift_c, time_reversal = \
            self.construct_symmetry_operators(K)
        A_cv = wfs.gd.cell_cv
        M_vv = np.dot(np.dot(A_cv.T, U_cc.T), np.linalg.inv(A_cv).T)
        ik = wfs.kd.bz2ibz_k[K]
        kpt = wfs.kpt_u[s * wfs.kd.nibzkpts + ik]
        psit_nG = kpt.psit_nG
        iG_Gv = 1j * wfs.pd.get_reciprocal_vectors(q=ik, add_q=False)
        ut_nvR = wfs.gd.zeros((n2 - n1, 3), complex)
        for n in range(n1, n2):
            for v in range(3):
                if self.real_space_derivatives:
                    ut_R = T(wfs.pd.ifft(psit_nG[n], ik))
                    grad_v[v](ut_R, ut_nvR[n - n1, v],
                              np.ones((3, 2), complex))
                else:
                    ut_R = T(wfs.pd.ifft(iG_Gv[:, v] * psit_nG[n], ik))
                    for v2 in range(3):
                        ut_nvR[n - n1, v2] += ut_R * M_vv[v, v2]

        return ut_nvR
