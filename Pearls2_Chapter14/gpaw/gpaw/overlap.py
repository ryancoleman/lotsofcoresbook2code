# -*- coding: utf-8 -*-
# Copyright (C) 2008  CSC Scientific Computing Ltd.
# Please see the accompanying LICENSE file for further information.

"""This module defines an Overlap operator.

The module defines an overlap operator and implements overlap-related
functions.

"""
import sys
import numpy as np

from gpaw.utilities.tools import lowdin, tri2full
from gpaw import extra_parameters
from gpaw.utilities.lapack import diagonalize
from gpaw import use_mic
from gpaw.mic import stream

class Overlap:
    """Overlap operator S

    This class contains information required to apply the
    overlap operator to a set of wavefunctions.
    """

    def __init__(self, ksl, timer):
        """Create the Overlap operator."""
        self.ksl = ksl
        self.timer = timer
        
    def orthonormalize(self, wfs, kpt, psit_nG=None):
        """Orthonormalizes the vectors a_nG with respect to the overlap.

        First, a Cholesky factorization C is done for the overlap
        matrix S_nn = <a_nG | S | a_nG> = C*_nn C_nn Cholesky matrix C
        is inverted and orthonormal vectors a_nG' are obtained as::

          psit_nG' = inv(C_nn) psit_nG
                    __
           ~   _   \    -1   ~   _
          psi (r) = )  C    psi (r)
             n     /__  nm     m
                    m

        Parameters
        ----------

        psit_nG: ndarray, input/output
            On input the set of vectors to orthonormalize,
            on output the overlap-orthonormalized vectors.
        kpt: KPoint object:
            k-point object from kpoint.py.
        work_nG: ndarray
            Optional work array for overlap matrix times psit_nG.
        work_nn: ndarray
            Optional work array for overlap matrix.

        """
        self.timer.start('Orthonormalize')
        if psit_nG is None:
            psit_nG = kpt.psit_nG
            if use_mic:
                psit_nG_mic = kpt.psit_nG_mic
        else:
            if use_mic:
                psit_nG_mic = stream.bind(psit_nG, update_device=False)
                stream.sync()

        P_ani = kpt.P_ani
        self.timer.start('projections')
        wfs.pt.integrate(psit_nG, P_ani, kpt.q)
        self.timer.stop('projections')

        # Construct the overlap matrix:
        operator = wfs.matrixoperator

        def S(psit_G):
            return psit_G
        
        def dS(a, P_ni):
            return np.dot(P_ni, wfs.setups[a].dO_ii)

        if use_mic:
            self.timer.start('calc_s_matrix')
            psit_nG_mic.update_device()
            stream.sync()
            S_nn = operator.calculate_matrix_elements(psit_nG_mic, P_ani, S, dS)
            self.timer.stop('calc_s_matrix')
        else:
            self.timer.start('calc_s_matrix')
            S_nn = operator.calculate_matrix_elements(psit_nG, P_ani, S, dS)
            self.timer.stop('calc_s_matrix')


        orthonormalization_string = repr(self.ksl)
        self.timer.start(orthonormalization_string)
        #
        if extra_parameters.get('sic', False):
            #
            # symmetric Loewdin Orthonormalization
            tri2full(S_nn, UL='L', map=np.conj)
            nrm_n = np.empty(S_nn.shape[0])
            diagonalize(S_nn, nrm_n)
            nrm_nn = np.diag(1.0/np.sqrt(nrm_n))
            S_nn = np.dot(np.dot(S_nn.T.conj(), nrm_nn), S_nn)
        else:
            #
            self.ksl.inverse_cholesky(S_nn)
        # S_nn now contains the inverse of the Cholesky factorization.
        # Let's call it something different:
        C_nn = S_nn
        del S_nn
        self.timer.stop(orthonormalization_string)

        self.timer.start('rotate_psi')
        if use_mic:
            operator.matrix_multiply(C_nn, psit_nG_mic, P_ani, out_nG=kpt.psit_nG_mic)
            kpt.psit_nG_mic.update_host()
            stream.sync()
            # kpt.psit_nG[:] = self.psit_nG_mic.array[:]
        else:
            operator.matrix_multiply(C_nn, psit_nG, P_ani, out_nG=kpt.psit_nG)
        self.timer.stop('rotate_psi')
        self.timer.stop('Orthonormalize')

    def apply(self, a_xG, b_xG, wfs, kpt, calculate_P_ani=True):
        """Apply the overlap operator to a set of vectors.

        Parameters
        ==========
        a_nG: ndarray
            Set of vectors to which the overlap operator is applied.
        b_nG: ndarray, output
            Resulting S times a_nG vectors.
        kpt: KPoint object
            k-point object defined in kpoint.py.
        calculate_P_ani: bool
            When True, the integrals of projector times vectors
            P_ni = <p_i | a_nG> are calculated.
            When False, existing P_ani are used

        """
        self.timer.start('Apply overlap')
        b_xG[:] = a_xG
        shape = a_xG.shape[:-3]
        P_axi = wfs.pt.dict(shape)

        if calculate_P_ani:
            wfs.pt.integrate(a_xG, P_axi, kpt.q)
        else:
            for a, P_ni in kpt.P_ani.items():
                P_axi[a][:] = P_ni

        for a, P_xi in P_axi.items():
            P_axi[a] = np.dot(P_xi, wfs.setups[a].dO_ii)
            # gemm(1.0, wfs.setups[a].dO_ii, P_xi, 0.0, P_xi, 'n')
        wfs.pt.add(b_xG, P_axi, kpt.q) # b_xG += sum_ai pt^a_i P_axi
        self.timer.stop('Apply overlap')

    def apply_inverse(self, a_xG, b_xG, wfs, kpt, calculate_P_ani=True):
        """Apply approximative inverse overlap operator to wave functions."""

        b_xG[:] = a_xG
        shape = a_xG.shape[:-3]
        P_axi = wfs.pt.dict(shape)

        if calculate_P_ani:
            wfs.pt.integrate(a_xG, P_axi, kpt.q)
        else:
            for a,P_ni in kpt.P_ani.items():
                P_axi[a][:] = P_ni

        for a, P_xi in P_axi.items():
            P_axi[a] = np.dot(P_xi, wfs.setups[a].dC_ii)
        wfs.pt.add(b_xG, P_axi, kpt.q)

