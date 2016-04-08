"""Module defining  ``Eigensolver`` classes."""

import numpy as np

from gpaw.utilities.blas import axpy
from gpaw.eigensolvers.eigensolver import Eigensolver

from gpaw import extra_parameters

class RMM_DIIS(Eigensolver):
    """RMM-DIIS eigensolver

    It is expected that the trial wave functions are orthonormal
    and the integrals of projector functions and wave functions
    ``nucleus.P_uni`` are already calculated

    Solution steps are:

    * Subspace diagonalization
    * Calculation of residuals
    * Improvement of wave functions:  psi' = psi + lambda PR + lambda PR'
    * Orthonormalization"""

    def __init__(self, keep_htpsit=True, blocksize=10,
                 fixed_trial_step=None):
        self.fixed_trial_step = fixed_trial_step
        Eigensolver.__init__(self, keep_htpsit, blocksize)

    def iterate_one_k_point(self, hamiltonian, wfs, kpt):
        """Do a single RMM-DIIS iteration for the kpoint"""

        psit_nG, R_nG = self.subspace_diagonalize(hamiltonian, wfs, kpt)

        self.timer.start('RMM-DIIS')
        if self.keep_htpsit:
            self.calculate_residuals(kpt, wfs, hamiltonian, psit_nG,
                                     kpt.P_ani, kpt.eps_n, R_nG)

        def integrate(a_G, b_G):
            return np.real(wfs.integrate(a_G, b_G, global_integral=False))

        comm = wfs.gd.comm
        B = self.blocksize
        dR_xG = wfs.empty(B, q=kpt.q)
        P_axi = wfs.pt.dict(B)
        error = 0.0
        for n1 in range(0, wfs.bd.mynbands, B):
            n2 = n1 + B
            if n2 > wfs.bd.mynbands:
                n2 = wfs.bd.mynbands
                B = n2 - n1
                P_axi = dict((a, P_xi[:B]) for a, P_xi in P_axi.items())
                dR_xG = dR_xG[:B]
                
            n_x = range(n1, n2)
            psit_xG = psit_nG[n1:n2]
            
            if self.keep_htpsit:
                R_xG = R_nG[n1:n2]
            else:
                R_xG = wfs.empty(B, q=kpt.q)
                wfs.apply_pseudo_hamiltonian(kpt, hamiltonian, psit_xG, R_xG)
                wfs.pt.integrate(psit_xG, P_axi, kpt.q)
                self.calculate_residuals(kpt, wfs, hamiltonian, psit_xG,
                                         P_axi, kpt.eps_n[n_x], R_xG, n_x)

            for n in n_x:
                if kpt.f_n is None:
                    weight = kpt.weight
                else:
                    weight = kpt.f_n[n]
                if self.nbands_converge != 'occupied':
                    if wfs.bd.global_index(n) < self.nbands_converge:
                        weight = kpt.weight
                    else:
                        weight = 0.0
                error += weight * integrate(R_xG[n - n1], R_xG[n - n1])

            # Precondition the residual:
            self.timer.start('precondition')
            ekin_x = self.preconditioner.calculate_kinetic_energy(
                psit_xG, kpt)
            dpsit_xG = self.preconditioner(R_xG, kpt, ekin_x)
            self.timer.stop('precondition')

            # Calculate the residual of dpsit_G, dR_G = (H - e S) dpsit_G:
            wfs.apply_pseudo_hamiltonian(kpt, hamiltonian, dpsit_xG, dR_xG)
            self.timer.start('projections')
            wfs.pt.integrate(dpsit_xG, P_axi, kpt.q)
            self.timer.stop('projections')
            self.calculate_residuals(kpt, wfs, hamiltonian, dpsit_xG,
                                     P_axi, kpt.eps_n[n_x], dR_xG, n_x,
                                     calculate_change=True)
            
            # Find lam that minimizes the norm of R'_G = R_G + lam dR_G
            RdR_x = np.array([integrate(dR_G, R_G)
                              for R_G, dR_G in zip(R_xG, dR_xG)])
            dRdR_x = np.array([integrate(dR_G, dR_G) for dR_G in dR_xG])
            comm.sum(RdR_x)
            comm.sum(dRdR_x)

            lam_x = -RdR_x / dRdR_x
            if extra_parameters.get('PK', False):
                lam_x[:] = np.where(lam_x>0.0, lam_x, 0.2)   
            # Calculate new psi'_G = psi_G + lam pR_G + lam2 pR'_G
            #                      = psi_G + p((lam+lam2) R_G + lam*lam2 dR_G)
            for lam, R_G, dR_G in zip(lam_x, R_xG, dR_xG):
                if self.fixed_trial_step is None:
                    lam2 = lam
                else:
                    lam2 = self.fixed_trial_step
                R_G *= lam + lam2
                axpy(lam * lam2, dR_G, R_G)
                
            self.timer.start('precondition')
            psit_xG[:] += self.preconditioner(R_xG, kpt, ekin_x)
            self.timer.stop('precondition')
            
        self.timer.stop('RMM-DIIS')
        error = comm.sum(error)
        return error, psit_nG
