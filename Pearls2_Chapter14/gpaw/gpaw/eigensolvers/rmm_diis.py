"""Module defining  ``Eigensolver`` classes."""

import numpy as np

from gpaw.utilities.blas import axpy, gemm, rk, r2k
from gpaw.utilities.lapack import general_diagonalize
from gpaw.utilities import unpack
from gpaw.eigensolvers.eigensolver import Eigensolver


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

    def __init__(self, keep_htpsit=True, blocksize=10, niter=3, rtol=1e-16,
                 limit_lambda=False, use_rayleigh=False, trial_step=0.1):
        """Initialize RMM-DIIS eigensolver.

        Parameters:

        limit_lambda: dictionary
            determines if step length should be limited
            supported keys:
            'absolute':True/False limit the absolute value
            'upper':float upper limit for lambda
            'lower':float lower limit for lambda

        """

        Eigensolver.__init__(self, keep_htpsit, blocksize)
        self.niter = niter
        self.rtol = rtol
        self.limit_lambda = limit_lambda
        self.use_rayleigh = use_rayleigh
        if use_rayleigh:
            self.blocksize = 1
        self.trial_step = trial_step
        self.first = True

    def iterate_one_k_point(self, hamiltonian, wfs, kpt):
        """Do a single RMM-DIIS iteration for the kpoint"""

        psit_nG, R_nG = self.subspace_diagonalize(hamiltonian, wfs, kpt)

        self.timer.start('RMM-DIIS')
        self.timer.start('Calculate residuals')
        if self.keep_htpsit:
            self.calculate_residuals(kpt, wfs, hamiltonian, psit_nG,
                                     kpt.P_ani, kpt.eps_n, R_nG)
        self.timer.stop('Calculate residuals')

        def integrate(a_G, b_G):
            return np.real(wfs.integrate(a_G, b_G, global_integral=False))

        comm = wfs.gd.comm
        B = self.blocksize
        dR_xG = wfs.empty(B, q=kpt.q)
        P_axi = wfs.pt.dict(B)
        errors_x = np.zeros(B)
        state_done = np.zeros(B, dtype=bool)

        errors_n = np.zeros(wfs.bd.mynbands)
        # Arrays needed for DIIS step
        if self.niter > 1:
            psit_diis_nxG = wfs.empty(B * self.niter, q=kpt.q)
            R_diis_nxG = wfs.empty(B * self.niter, q=kpt.q)
            # P_diis_anxi = wfs.pt.dict(B * self.niter)
            eig_n = np.zeros(self.niter)  # eigenvalues for diagonalization
                                          # not needed in any step

        error = 0.0
        for n1 in range(0, wfs.bd.mynbands, B):
            state_done[:] = False
            n2 = n1 + B
            if n2 > wfs.bd.mynbands:
                n2 = wfs.bd.mynbands
                B = n2 - n1
                P_axi = dict((a, P_xi[:B]) for a, P_xi in P_axi.items())
                dR_xG = dR_xG[:B]
                
            n_x = range(n1, n2)
            psit_xG = psit_nG[n1:n2]

            self.timer.start('Calculate residuals')
            if self.keep_htpsit:
                R_xG = R_nG[n1:n2]
            else:
                R_xG = wfs.empty(B, q=kpt.q)
                wfs.apply_pseudo_hamiltonian(kpt, hamiltonian, psit_xG, R_xG)
                wfs.pt.integrate(psit_xG, P_axi, kpt.q)
                self.calculate_residuals(kpt, wfs, hamiltonian, psit_xG,
                                         P_axi, kpt.eps_n[n_x], R_xG, n_x)
            self.timer.stop('Calculate residuals')

            errors_x[:] = 0.0
            for n in range(n1, n2):
                if kpt.f_n is None:
                    weight = kpt.weight
                else:
                    weight = kpt.f_n[n]
                if self.nbands_converge != 'occupied':
                    if wfs.bd.global_index(n) < self.nbands_converge:
                        weight = kpt.weight
                    else:
                        weight = 0.0
                errors_x[n - n1] = weight * integrate(R_xG[n - n1], 
                                                  R_xG[n - n1])
                errors_n[n] = errors_x[n - n1]
            comm.sum(errors_x)
            error += np.sum(errors_x)

            # Insert first vectors and residuals for DIIS step
            if self.niter > 1:
                # Save the previous vectors contiguously for each band
                # in the block
                psit_diis_nxG[:B * self.niter:self.niter] = psit_xG
                R_diis_nxG[:B * self.niter:self.niter] = R_xG

            # Precondition the residual:
            self.timer.start('precondition')
            # ekin_x = self.preconditioner.calculate_kinetic_energy(
            #     R_xG, kpt)
            ekin_x = self.preconditioner.calculate_kinetic_energy(
                psit_xG, kpt)
            dpsit_xG = self.preconditioner(R_xG, kpt, ekin_x)
            self.timer.stop('precondition')

            # Calculate the residual of dpsit_G, dR_G = (H - e S) dpsit_G:
            # self.timer.start('Apply Hamiltonian')
            wfs.apply_pseudo_hamiltonian(kpt, hamiltonian, dpsit_xG, dR_xG)
            # self.timer.stop('Apply Hamiltonian')
            self.timer.start('projections')
            wfs.pt.integrate(dpsit_xG, P_axi, kpt.q)
            self.timer.stop('projections')

            if self.use_rayleigh:
                self.timer.start('Minimize Rayleigh')
                i1 = wfs.integrate(psit_xG, dR_xG, 
                                          global_integral=False).item()
                i2 = wfs.integrate(dpsit_xG, dR_xG, 
                                          global_integral=False).item()
                i3 = wfs.integrate(dpsit_xG, psit_xG, 
                                          global_integral=False).item()
                i4 = wfs.integrate(dpsit_xG, dpsit_xG, 
                                          global_integral=False).item()

                for a, dP_xi in P_axi.items():
                    P_i = kpt.P_ani[a][n1]
                    dP_i = dP_xi[0]
                    dH_ii = unpack(hamiltonian.dH_asp[a][kpt.s])
                    dO_ii = wfs.setups[a].dO_ii
                    i1 += np.dot(P_i, 
                                        np.dot(dH_ii, dP_i.conj())).item()
                    i2 += np.dot(dP_i, 
                                        np.dot(dH_ii, dP_i.conj())).item()
                    i3 += np.dot(dP_i, 
                                        np.dot(dO_ii, P_i.conj())).item()
                    i4 += np.dot(dP_i, 
                                        np.dot(dO_ii, dP_i.conj())).item()

                i1 = comm.sum(i1)
                i2 = comm.sum(i2)
                i3 = comm.sum(i3)
                i4 = comm.sum(i4)

                a = np.real(i2 * i3 - i1 * i4)
                b = np.real(i2 - kpt.eps_n[n1] * i4)
                c = np.real(i1 - kpt.eps_n[n1] * i3)
                # print "A,B,C", a,b,c
                lam_x = np.array((-2.0 * c / (b + np.sqrt(b**2 - 4.0 * a * c)),))

                self.timer.stop('Minimize Rayleigh')
                self.timer.start('Calculate residuals')
                self.calculate_residuals(kpt, wfs, hamiltonian, dpsit_xG,
                                         P_axi, kpt.eps_n[n_x], dR_xG, n_x,
                                         calculate_change=True)
                self.timer.stop('Calculate residuals')
            else:
                self.timer.start('Calculate residuals')
                self.calculate_residuals(kpt, wfs, hamiltonian, dpsit_xG,
                                         P_axi, kpt.eps_n[n_x], dR_xG, n_x,
                                         calculate_change=True)
                self.timer.stop('Calculate residuals')

                # Find lam that minimizes the norm of R'_G = R_G + lam dR_G
                self.timer.start('Find lambda')
                RdR_x = np.array([integrate(dR_G, R_G)
                                  for R_G, dR_G in zip(R_xG, dR_xG)])
                dRdR_x = np.array([integrate(dR_G, dR_G) for dR_G in dR_xG])
                comm.sum(RdR_x)
                comm.sum(dRdR_x)
                lam_x = -RdR_x / dRdR_x
                self.timer.stop('Find lambda')

            # print "Lam_x:", lam_x
            # Limit abs(lam) to [0.15, 1.0]
            if self.limit_lambda:
                upper = self.limit_lambda['upper']
                lower = self.limit_lambda['lower']
                if self.limit_lambda.get('absolute', False):
                    lam_x = np.where(np.abs(lam_x) < lower, 
                                     lower * np.sign(lam_x), lam_x)
                    lam_x = np.where(np.abs(lam_x) > upper, 
                                     upper * np.sign(lam_x), lam_x)
                else:
                    lam_x = np.where(lam_x < lower, lower, lam_x)
                    lam_x = np.where(lam_x > upper, upper, lam_x)

            # lam_x[:] = 0.1

            # New trial wavefunction and residual          
            self.timer.start('Update psi')
            for lam, psit_G, dpsit_G, R_G, dR_G in zip(lam_x, psit_xG, 
                                                       dpsit_xG, R_xG, dR_xG):
                axpy(lam, dpsit_G, psit_G)  # psit_G += lam * dpsit_G
                axpy(lam, dR_G, R_G)  # R_G += lam** dR_G
            self.timer.stop('Update psi')

            self.timer.start('DIIS step')
            # DIIS step
            for nit in range(1, self.niter):
                # Do not perform DIIS if error is small
                # if abs(error_block / B) < self.rtol:
                #     break 
                
                # Update the subspace
                psit_diis_nxG[nit:B * self.niter:self.niter] = psit_xG
                R_diis_nxG[nit:B * self.niter:self.niter] = R_xG

                # XXX Only integrals of nit old psits would be needed
                # self.timer.start('projections')
                # wfs.pt.integrate(psit_diis_nxG, P_diis_anxi, kpt.q)
                # self.timer.stop('projections')
                if nit > 1 or self.limit_lambda:
                    for ib in range(B):
                        if state_done[ib]:
                            continue
                        istart = ib * self.niter
                        iend = istart + nit + 1

                        # Residual matrix
                        self.timer.start('Construct matrix')
                        R_nn = wfs.integrate(R_diis_nxG[istart:iend],
                                             R_diis_nxG[istart:iend],
                                             global_integral=True)

                        # Full matrix                     
                        A_nn = -np.ones((nit + 2, nit + 2), wfs.dtype)
                        A_nn[:nit+1, :nit+1] = R_nn[:]
                        A_nn[-1,-1] = 0.0
                        x_n = np.zeros(nit + 2, wfs.dtype)
                        x_n[-1] = -1.0
                        self.timer.stop('Construct matrix')
                        self.timer.start('Linear solve')
                        alpha_i = np.linalg.solve(A_nn, x_n)[:-1]
                        self.timer.stop('Linear solve')
                        self.timer.start('Update trial vectors')
                        psit_xG[ib] = alpha_i[nit] * psit_diis_nxG[istart + nit]
                        R_xG[ib] = alpha_i[nit] * R_diis_nxG[istart + nit]
                        for i in range(nit):
                            # axpy(alpha_i[i], psit_diis_nxG[istart + i], 
                            #      psit_diis_nxG[istart + nit])
                            # axpy(alpha_i[i], R_diis_nxG[istart + i], 
                            #      R_diis_nxG[istart + nit])
                            axpy(alpha_i[i], psit_diis_nxG[istart + i], 
                                 psit_xG[ib])
                            axpy(alpha_i[i], R_diis_nxG[istart + i], 
                                 R_xG[ib])
                        self.timer.stop('Update trial vectors')

                if nit < self.niter - 1:
                    self.timer.start('precondition')
                    # ekin_x = self.preconditioner.calculate_kinetic_energy(
                    #     R_xG, kpt)
                    dpsit_xG = self.preconditioner(R_xG, kpt, ekin_x)
                    self.timer.stop('precondition')

                    for psit_G, lam, dpsit_G in zip(psit_xG, lam_x, dpsit_xG):
                        axpy(lam, dpsit_G, psit_G)

                    # Calculate the new residuals
                    self.timer.start('Calculate residuals')
                    wfs.apply_pseudo_hamiltonian(kpt, hamiltonian, psit_xG, 
                                                 R_xG)
                    wfs.pt.integrate(psit_xG, P_axi, kpt.q)
                    self.calculate_residuals(kpt, wfs, hamiltonian, psit_xG,
                                             P_axi, kpt.eps_n[n_x], R_xG, n_x,
                                             calculate_change=True)
                    self.timer.stop('Calculate residuals')
                    self.timer.start('Calculate errors')                
                    errors_new_x = np.zeros(B)
                    # errors_x[:] = 0.0
                    for n in range(n1, n2):
                        if kpt.f_n is None:
                            weight = kpt.weight
                        else:
                            weight = kpt.f_n[n]
                        if self.nbands_converge != 'occupied':
                            if wfs.bd.global_index(n) < self.nbands_converge:
                                weight = kpt.weight
                            else:
                                weight = 0.0
                        errors_new_x[n-n1] += weight * integrate(R_xG[n - n1], 
                                                      R_xG[n - n1])
                    comm.sum(errors_x)
                    self.timer.stop('Calculate errors')                

            self.timer.stop('DIIS step')                
            # Final trial step
            self.timer.start('precondition')
            # ekin_x = self.preconditioner.calculate_kinetic_energy(
            #     R_xG, kpt)
            dpsit_xG = self.preconditioner(R_xG, kpt, ekin_x)
            self.timer.stop('precondition')
            self.timer.start('Update psi')
                
            if self.trial_step is not None:
                lam_x[:] = self.trial_step
            for lam, psit_G, dpsit_G in zip(lam_x, psit_xG, dpsit_xG):
                axpy(lam, dpsit_G, psit_G)  # psit_G += lam * dpsit_G
            self.timer.stop('Update psi')

            # norm = wfs.integrate(psit_xG[0], psit_xG[0])
            # wfs.pt.integrate(psit_xG, P_axi, kpt.q)
            # for a, P_xi in P_axi.items():
            #     dO_ii = wfs.setups[a].dO_ii
            #     norm += np.vdot(P_xi[0], np.inner(dO_ii, P_xi[0]))
            # norm = comm.sum(np.real(norm).item())
            # psit_xG /= np.sqrt(norm)
            
        self.timer.stop('RMM-DIIS')
        return error, psit_nG

    def __repr__(self):
        repr_string = 'RMM-DIIS eigensolver\n'
        repr_string += '       keep_htpsit: %s\n' % self.keep_htpsit
        repr_string += '       Block size: %d\n' % self.blocksize
        repr_string += '       DIIS iterations: %d\n' % self.niter
        repr_string += '       Threshold for DIIS: %5.1e\n' % self.rtol
        repr_string += '       Limit lambda: %s\n' % self.limit_lambda
        repr_string += '       use_rayleigh: %s\n' % self.use_rayleigh
        repr_string += '       trial_step: %s' % self.trial_step
        return repr_string

