"""Module defining an eigensolver base-class."""

import numpy as np

from gpaw.utilities.blas import axpy
from gpaw.utilities import unpack
from gpaw.hs_operators import reshape
from gpaw import use_mic
from gpaw.mic import stream


class Eigensolver:
    def __init__(self, keep_htpsit=True, blocksize=1):
        self.keep_htpsit = keep_htpsit
        self.initialized = False
        self.Htpsit_nG = None
        self.error = np.inf
        self.blocksize = blocksize
        self.orthonormalization_required = True
        
    def initialize(self, wfs):
        self.timer = wfs.timer
        self.world = wfs.world
        self.kpt_comm = wfs.kd.comm
        self.band_comm = wfs.band_comm
        self.dtype = wfs.dtype
        self.bd = wfs.bd
        self.ksl = wfs.diagksl
        self.nbands = wfs.bd.nbands
        self.mynbands = wfs.bd.mynbands
        self.operator = wfs.matrixoperator

        if self.mynbands != self.nbands or self.operator.nblocks != 1:
            self.keep_htpsit = False

        if self.keep_htpsit:
            self.Htpsit_nG = wfs.empty(self.nbands)
            if use_mic:
                self.Htpsit_nG_mic = stream.bind(self.Htpsit_nG)
                stream.sync()

        # Preconditioner for the electronic gradients:
        self.preconditioner = wfs.make_preconditioner(self.blocksize)

        for kpt in wfs.kpt_u:
            if kpt.eps_n is None:
                kpt.eps_n = np.empty(self.mynbands)
        
        # Allocate arrays for matrix operator
        self.operator.allocate_arrays()

        self.initialized = True

    def reset(self):
        self.initialized = False

    def iterate(self, hamiltonian, wfs):
        """Solves eigenvalue problem iteratively

        This method is inherited by the actual eigensolver which should
        implement *iterate_one_k_point* method for a single iteration of
        a single kpoint.
        """

        if not self.initialized:
            self.initialize(wfs)

        error = 0.0
        for kpt in wfs.kpt_u:
            if not wfs.orthonormalized:
                wfs.overlap.orthonormalize(wfs, kpt)
            e, psit_nG = self.iterate_one_k_point(hamiltonian, wfs, kpt)
            error += e
            if self.orthonormalization_required:
                wfs.overlap.orthonormalize(wfs, kpt, psit_nG)
            else:
                kpt.psit_nG[:] = psit_nG[:]

        wfs.set_orthonormalized(True)

        self.error = self.band_comm.sum(self.kpt_comm.sum(error))

    def iterate_one_k_point(self, hamiltonian, kpt):
        """Implemented in subclasses."""
        raise NotImplementedError

    def calculate_residuals(self, kpt, wfs, hamiltonian, psit_xG, P_axi, eps_x,
                            R_xG, n_x=None, calculate_change=False):
        """Calculate residual.

        From R=Ht*psit calculate R=H*psit-eps*S*psit."""
        
        for R_G, eps, psit_G in zip(R_xG, eps_x, psit_xG):
            axpy(-eps, psit_G, R_G)

        c_axi = {}
        for a, P_xi in P_axi.items():
            dH_ii = unpack(hamiltonian.dH_asp[a][kpt.s])
            dO_ii = hamiltonian.setups[a].dO_ii
            c_xi = (np.dot(P_xi, dH_ii) -
                    np.dot(P_xi * eps_x[:, np.newaxis], dO_ii))
            c_axi[a] = c_xi
        hamiltonian.xc.add_correction(kpt, psit_xG, R_xG, P_axi, c_axi, n_x,
                                      calculate_change)
        wfs.pt.add(R_xG, c_axi, kpt.q)
        
    def subspace_diagonalize(self, hamiltonian, wfs, kpt):
        """Diagonalize the Hamiltonian in the subspace of kpt.psit_nG

        *Htpsit_nG* is a work array of same size as psit_nG which contains
        the local part of the Hamiltonian times psit on exit

        First, the Hamiltonian (defined by *kin*, *vt_sG*, and
        *dH_asp*) is applied to the wave functions, then the *H_nn*
        matrix is calculated and diagonalized, and finally, the wave
        functions (and also Htpsit_nG are rotated.  Also the
        projections *P_ani* are rotated.

        It is assumed that the wave functions *psit_nG* are orthonormal
        and that the integrals of projector functions and wave functions
        *P_ani* are already calculated.

        Return ratated wave functions and H applied to the rotated
        wave functions if self.keep_htpsit is True.
        """

        if self.band_comm.size > 1 and wfs.bd.strided:
            raise NotImplementedError

        self.timer.start('Subspace diag')

        if use_mic:
            psit_nG = kpt.psit_nG_mic
            # psit_nG.update_device()
            # stream.sync()
        else:
            psit_nG = kpt.psit_nG
        P_ani = kpt.P_ani

        if self.keep_htpsit:
            if use_mic:
                Htpsit_nG = self.Htpsit_nG_mic
            else:
                Htpsit_nG = reshape(self.Htpsit_nG, psit_nG.shape)
        else:
            Htpsit_nG = None

        def H(psit_xG):
            if self.keep_htpsit:
                result_xG = Htpsit_nG
            else:
                if use_mic:
                    result_xG = self.operator.work1_xG_mic
                else:
                    result_xG = reshape(self.operator.work1_xG, psit_xG.shape)
            if use_mic:
                psit_xG.update_device()
                wfs.apply_pseudo_hamiltonian(kpt, hamiltonian, psit_xG.array,
                                             result_xG.array)
                result_xG.update_device() 
                stream.sync()
            else:
                wfs.apply_pseudo_hamiltonian(kpt, hamiltonian, psit_xG,
                                             result_xG)
            hamiltonian.xc.apply_orbital_dependent_hamiltonian(
                kpt, psit_xG, result_xG, hamiltonian.dH_asp)
            return result_xG

        def dH(a, P_ni):
            return np.dot(P_ni, unpack(hamiltonian.dH_asp[a][kpt.s]))

        self.timer.start('calc_h_matrix')
        H_nn = self.operator.calculate_matrix_elements(psit_nG, P_ani,
                                                       H, dH)
        hamiltonian.xc.correct_hamiltonian_matrix(kpt, H_nn)
        self.timer.stop('calc_h_matrix')

        diagonalization_string = repr(self.ksl)
        wfs.timer.start(diagonalization_string)
        self.ksl.diagonalize(H_nn, kpt.eps_n)
        # H_nn now contains the result of the diagonalization.
        wfs.timer.stop(diagonalization_string)

        self.timer.start('rotate_psi')
        psit_nG = self.operator.matrix_multiply(H_nn, psit_nG, P_ani)
        if self.keep_htpsit:
            if use_mic:
                Htpsit_nG = self.operator.matrix_multiply(H_nn, Htpsit_nG,
                                                          out_nG=kpt.psit_nG_mic)
                 
            else:
                Htpsit_nG = self.operator.matrix_multiply(H_nn, Htpsit_nG,
                                                          out_nG=kpt.psit_nG)

        # Rotate orbital dependent XC stuff:
        hamiltonian.xc.rotate(kpt, H_nn)

        self.timer.stop('rotate_psi')
        self.timer.stop('Subspace diag')

        if use_mic:
            psit_nG.update_host()
            stream.sync()
            if self.keep_htpsit:
                Htpsit_nG.update_host()
                stream.sync()
                return psit_nG.array, Htpsit_nG.array
            else:
                return psit_nG.array, Htpsit_nG
        else:
            return psit_nG, Htpsit_nG

    def estimate_memory(self, mem, wfs):
        gridmem = wfs.bytes_per_wave_function()

        keep_htpsit = self.keep_htpsit and (wfs.bd.mynbands == wfs.bd.nbands)

        if keep_htpsit:
            mem.subnode('Htpsit', wfs.bd.nbands * gridmem)
        else:
            mem.subnode('No Htpsit', 0)

        mem.subnode('eps_n', wfs.bd.mynbands * mem.floatsize)
        mem.subnode('eps_N', wfs.bd.nbands * mem.floatsize)
        mem.subnode('Preconditioner', 4 * gridmem)
        mem.subnode('Work', gridmem)
