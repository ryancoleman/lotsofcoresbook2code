"""Module defining  ``Eigensolver`` classes."""

from math import pi, sqrt, sin, cos, atan2

import numpy as np
from numpy import dot
from ase.units import Hartree

from gpaw.utilities.blas import axpy, dotc, gemv, gemm
from gpaw.utilities import unpack
from gpaw.eigensolvers.eigensolver import Eigensolver
from gpaw.hs_operators import reshape

from gpaw import extra_parameters

class CG(Eigensolver):
    """Conjugate gardient eigensolver

    It is expected that the trial wave functions are orthonormal
    and the integrals of projector functions and wave functions
    ``nucleus.P_uni`` are already calculated.

    Solution steps are:

    * Subspace diagonalization
    * Calculate all residuals
    * Conjugate gradient steps
    """

    def __init__(self, niter=4, rtol=0.30):
        """Construct conjugate gradient eigen solver.

        parameters:

        niter: int
            Maximum number of conjugate gradient iterations per band
        rtol: float
            If change in residual is less than rtol, iteration for band is
            not continued

        """
        Eigensolver.__init__(self)
        self.niter = niter
        self.rtol = rtol
        if extra_parameters.get('PK', False):
            self.orthonormalization_required = True
        else:
            self.orthonormalization_required = False

    def __repr__(self):
        return 'CG(niter=%d, rtol=%5.1e)' % (self.niter, self.rtol)

    def initialize(self, wfs):
        if wfs.bd.comm.size > 1:
            raise ValueError('CG eigensolver does not support band '
                             'parallelization.  This calculation parallelizes '
                             'over %d band groups.' % wfs.bd.comm.size)
        Eigensolver.initialize(self, wfs)
        self.overlap = wfs.overlap

    def iterate_one_k_point(self, hamiltonian, wfs, kpt):
        """Do conjugate gradient iterations for the k-point"""

        niter = self.niter

        phi_G = wfs.empty(q=kpt.q)
        phi_old_G = wfs.empty(q=kpt.q)

        comm = wfs.gd.comm

        psit_nG, Htpsit_nG = self.subspace_diagonalize(hamiltonian, wfs, kpt)
        # Note that psit_nG is now in self.operator.work1_nG and
        # Htpsit_nG is in kpt.psit_nG!

        R_nG = reshape(self.Htpsit_nG, psit_nG.shape)
        Htphi_G = R_nG[0]

        R_nG[:] = Htpsit_nG
        self.timer.start('Residuals')
        self.calculate_residuals(kpt, wfs, hamiltonian, psit_nG,
                                 kpt.P_ani, kpt.eps_n, R_nG)
        self.timer.stop('Residuals')

        self.timer.start('CG')

        total_error = 0.0
        for n in range(self.nbands):
            if extra_parameters.get('PK', False):
               N = n+1
            else:
               N = psit_nG.shape[0]+1
            R_G = R_nG[n]
            Htpsit_G = Htpsit_nG[n]
            gamma_old = 1.0
            phi_old_G[:] = 0.0
            error = np.real(wfs.integrate(R_G, R_G))
            for nit in range(niter):
                if (error * Hartree**2 < self.tolerance / self.nbands):
                    break

                ekin = self.preconditioner.calculate_kinetic_energy(
                    psit_nG[n:n + 1], kpt)

                pR_G = self.preconditioner(R_nG[n:n + 1], kpt, ekin)

                # New search direction
                gamma = comm.sum(np.vdot(pR_G, R_G).real)
                phi_G[:] = -pR_G - gamma / gamma_old * phi_old_G
                gamma_old = gamma
                phi_old_G[:] = phi_G[:]

                # Calculate projections
                P2_ai = wfs.pt.dict()
                wfs.pt.integrate(phi_G, P2_ai, kpt.q)

                # Orthonormalize phi_G to all bands
                self.timer.start('CG: orthonormalize')
                self.timer.start('CG: overlap')
                overlap_n = wfs.integrate(psit_nG[:N], phi_G,
                                        global_integral=False)
                self.timer.stop('CG: overlap')
                self.timer.start('CG: overlap2')
                for a, P2_i in P2_ai.items():
                    P_ni = kpt.P_ani[a]
                    dO_ii = wfs.setups[a].dO_ii
                    gemv(1.0, P_ni[:N].conjugate(), np.inner(dO_ii, P2_i), 
                         1.0, overlap_n)
                self.timer.stop('CG: overlap2')
                comm.sum(overlap_n)

                # phi_G -= overlap_n * kpt.psit_nG
                wfs.matrixoperator.gd.gemv(-1.0, psit_nG[:N], overlap_n, 
                                           1.0, phi_G, 'n')
                for a, P2_i in P2_ai.items():
                    P_ni = kpt.P_ani[a]
                    gemv(-1.0, P_ni[:N], overlap_n, 1.0, P2_i, 'n')

                norm = wfs.integrate(phi_G, phi_G, global_integral=False)
                for a, P2_i in P2_ai.items():
                    dO_ii = wfs.setups[a].dO_ii
                    norm += np.vdot(P2_i, np.inner(dO_ii, P2_i))
                norm = comm.sum(np.real(norm).item())
                phi_G /= sqrt(norm)
                for P2_i in P2_ai.values():
                    P2_i /= sqrt(norm)
                self.timer.stop('CG: orthonormalize')

                # find optimum linear combination of psit_G and phi_G
                an = kpt.eps_n[n]
                wfs.apply_pseudo_hamiltonian(kpt, hamiltonian,
                                             phi_G.reshape((1,) + phi_G.shape),
                                             Htphi_G.reshape((1,) +
                                                             Htphi_G.shape))
                b = wfs.integrate(phi_G, Htpsit_G, global_integral=False)
                c = wfs.integrate(phi_G, Htphi_G, global_integral=False)
                for a, P2_i in P2_ai.items():
                    P_i = kpt.P_ani[a][n]
                    dH_ii = unpack(hamiltonian.dH_asp[a][kpt.s])
                    b += dot(P2_i, dot(dH_ii, P_i.conj()))
                    c += dot(P2_i, dot(dH_ii, P2_i.conj()))
                b = comm.sum(np.real(b).item())
                c = comm.sum(np.real(c).item())

                theta = 0.5 * atan2(2 * b, an - c)
                enew = (an * cos(theta)**2 +
                        c * sin(theta)**2 +
                        b * sin(2.0 * theta))
                # theta can correspond either minimum or maximum
                if (enew - kpt.eps_n[n]) > 0.0:  # we were at maximum
                    theta += pi / 2.0
                    enew = (an * cos(theta)**2 +
                            c * sin(theta)**2 +
                            b * sin(2.0 * theta))

                kpt.eps_n[n] = enew
                psit_nG[n] *= cos(theta)
                # kpt.psit_nG[n] += sin(theta) * phi_G
                axpy(sin(theta), phi_G, psit_nG[n])
                for a, P2_i in P2_ai.items():
                    P_i = kpt.P_ani[a][n]
                    P_i *= cos(theta)
                    P_i += sin(theta) * P2_i

                if nit < niter - 1:
                    Htpsit_G *= cos(theta)
                    # Htpsit_G += sin(theta) * Htphi_G
                    axpy(sin(theta), Htphi_G, Htpsit_G)
                    #adjust residuals
                    R_G[:] = Htpsit_G - kpt.eps_n[n] * psit_nG[n]

                    coef_ai = wfs.pt.dict()
                    for a, coef_i in coef_ai.items():
                        P_i = kpt.P_ani[a][n]
                        dO_ii = wfs.setups[a].dO_ii
                        dH_ii = unpack(hamiltonian.dH_asp[a][kpt.s])
                        coef_i[:] = (dot(P_i, dH_ii) -
                                     dot(P_i * kpt.eps_n[n], dO_ii))
                    wfs.pt.add(R_G, coef_ai, kpt.q)
                    error_new = np.real(wfs.integrate(R_G, R_G))
                    if error_new / error < self.rtol:
                        # print >> self.f, "cg:iters", n, nit+1
                        break
                    if (self.nbands_converge == 'occupied' and
                        kpt.f_n is not None and kpt.f_n[n] == 0.0):
                        # print >> self.f, "cg:iters", n, nit+1
                        break
                    error = error_new

            if kpt.f_n is None:
                weight = 1.0
            else:
                weight = kpt.f_n[n]
            if self.nbands_converge != 'occupied':
                weight = kpt.weight * float(n < self.nbands_converge)
            total_error += weight * error
            # if nit == 3:
            #   print >> self.f, "cg:iters", n, nit+1

        self.timer.stop('CG')
        return total_error, psit_nG
