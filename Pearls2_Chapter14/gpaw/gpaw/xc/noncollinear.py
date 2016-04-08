from math import sqrt, pi

import numpy as np

from gpaw.xc.functional import XCFunctional
from gpaw.xc.lda import LDA
from gpaw.xc.libxc import LibXC
from gpaw.lcao.eigensolver import DirectLCAO
from gpaw.wavefunctions.lcao import LCAOWaveFunctions
from gpaw.utilities import unpack
from gpaw.utilities.blas import gemm
from gpaw.mixer import BaseMixer
from gpaw.utilities.tools import tri2full
from gpaw.sphere.lebedev import Y_nL, weight_n
from gpaw.xc.pawcorrection import rnablaY_nLv

X = np.newaxis


class NonCollinearLDAKernel(LibXC):
    def __init__(self):
        LibXC.__init__(self, 'LDA')
        
    def calculate(self, e_g, n_sg, dedn_sg,
                  sigma_xg=None, dedsigma_xg=None,
                  tau_sg=None, dedtau_sg=None):
        n_g = n_sg[0]
        m_vg = n_sg[1:4]
        m_g = (m_vg**2).sum(0)**0.5
        nnew_sg = np.empty((2,) + n_g.shape)
        nnew_sg[:] = n_g
        nnew_sg[0] += m_g
        nnew_sg[1] -= m_g
        nnew_sg *= 0.5
        vnew_sg = np.zeros_like(nnew_sg)
        LibXC.calculate(self, e_g, nnew_sg, vnew_sg)
        dedn_sg[0] += 0.5 * vnew_sg.sum(0)
        vnew_sg /= np.where(m_g < 1e-15, 1, m_g)
        dedn_sg[1:4] += 0.5 * vnew_sg[0] * m_vg
        dedn_sg[1:4] -= 0.5 * vnew_sg[1] * m_vg


class NonCollinearFunctional(XCFunctional):
    def __init__(self, xc):
        XCFunctional.__init__(self, xc.name)
        self.xc = xc
        self.type = xc.type
        
    def calculate(self, gd, n_sg, dedn_sg=None, e_g=None):
        n_g = n_sg[0]
        m_vg = n_sg[1:4]
        m_g = (m_vg**2).sum(0)**0.5
        nnew_sg = gd.empty(2)
        nnew_sg[:] = n_g
        nnew_sg[0] += m_g
        nnew_sg[1] -= m_g
        nnew_sg *= 0.5
        vnew_sg = gd.zeros(2)
        if e_g is None:
            e_g = gd.empty()
        exc = self.xc.calculate(gd, nnew_sg, vnew_sg, e_g)
        if dedn_sg is not None:
            dedn_sg[0] += 0.5 * vnew_sg.sum(0)
            vnew_sg /= np.where(m_g < 1e-15, 1, m_g)
            dedn_sg[1:4] += 0.5 * vnew_sg[0] * m_vg
            dedn_sg[1:4] -= 0.5 * vnew_sg[1] * m_vg
        return exc
    
    def calculate_paw_correction(self, setup, D_sp, dEdD_sp=None,
                                 addcoredensity=True, a=None):
        c = setup.xc_correction
        if c is None:
            return 0.0
        
        assert addcoredensity
        assert len(D_sp) == 4
        assert c.nc_corehole_g is None
        
        rgd = c.rgd
        
        nc0_g = sqrt(4 * pi) * c.nc_g
        nct0_g = sqrt(4 * pi) * c.nct_g
        
        D_sLq = np.inner(D_sp, c.B_pqL.T)

        e, dEdD_sqL = self.calculate_radial_expansion(rgd, D_sLq, c.n_qg,
                                                      nc0_g)
        et, dEtdD_sqL = self.calculate_radial_expansion(rgd, D_sLq, c.nt_qg,
                                                        nct0_g)

        if dEdD_sp is not None:
            dEdD_sp += np.inner((dEdD_sqL - dEtdD_sqL).reshape((4, -1)),
                                c.B_pqL.reshape((len(c.B_pqL), -1)))

        return e - et - c.Exc0

    def calculate_radial_expansion(self, rgd, D_sLq, n_qg, nc0_g):
        D_Lq = D_sLq[0]
        M_vLq = D_sLq[1:]
        
        n_Lg = np.dot(D_Lq, n_qg)
        n_Lg[0] += nc0_g

        m_vLg = np.dot(M_vLq, n_qg)

        if self.xc.type == 'GGA':
            dndr_Lg = np.empty_like(n_Lg)
            rgd.derivative(n_Lg, dndr_Lg)
            dmdr_vLg = np.empty_like(m_vLg)
            rgd.derivative(m_vLg, dmdr_vLg)
        elif self.xc.type == 'LDA':
            dndr_Lg = None
            dmdr_vLg = None

        Lmax, nq = D_Lq.shape
        dEdD_sqL = np.zeros((4, nq, Lmax))
        
        E = 0.0
        for w, Y_L, rnablaY_Lv in zip(weight_n,
                                      Y_nL[:, :Lmax],
                                      rnablaY_nLv[:, :Lmax]):
            if self.xc.type == 'LDA':
                (e_g, dedn_g, dedm_vg) = \
                    self.calculate_radial(rgd, n_Lg, Y_L, dndr_Lg, rnablaY_Lv,
                                          m_vLg, dmdr_vLg)
            else:
                (e_g, dedn_g, dedm_vg,
                 a_g, b_vg, c_g, d_vg, dedsigma_xg, dcdm_vg, dddm_vLvg,
                 m_vg, m_g) = \
                 self.calculate_radial(rgd, n_Lg, Y_L, dndr_Lg, rnablaY_Lv,
                                       m_vLg, dmdr_vLg)
            dEdD_sqL[0] += np.dot(rgd.dv_g * dedn_g,
                                  n_qg.T)[:, X] * (w * Y_L)
            dEdD_sqL[1:] += np.dot(rgd.dv_g * dedm_vg,
                                   n_qg.T)[:, :, X] * (w * Y_L)
            if self.xc.type == 'GGA':
                v_g = rgd.empty()
                rgd.derivative2(rgd.dv_g * (dedsigma_xg[0] * (a_g + c_g) +
                                            dedsigma_xg[1] * a_g +
                                            dedsigma_xg[2] * (a_g - c_g)), v_g)
                dEdD_sqL[0] -= np.dot(v_g, n_qg.T)[:, X] * (0.5 * w * Y_L)
                B_vg = (dedsigma_xg[0] * (b_vg + d_vg) +
                        dedsigma_xg[1] * b_vg +
                        dedsigma_xg[2] * (b_vg - d_vg))
                B_vq = np.dot(B_vg * rgd.dr_g, n_qg.T)
                dEdD_sqL[0] += 2 * pi * w * np.dot(rnablaY_Lv, B_vq).T

                A_g = (dedsigma_xg[0] * (a_g + c_g) -
                       dedsigma_xg[1] * c_g +
                       dedsigma_xg[2] * (c_g - a_g)) * rgd.dv_g
                v_vg = rgd.empty(3)
                rgd.derivative2(A_g * m_vg / m_g, v_vg)
                v_vg -= A_g * dcdm_vg
                dEdD_sqL[1:] -= np.dot(v_vg, n_qg.T)[:, :, X] * (0.5 * w * Y_L)

                dedsigma_xg *= rgd.dr_g
                B_vg = (dedsigma_xg[0] * (b_vg + d_vg) -
                        dedsigma_xg[1] * d_vg +
                        dedsigma_xg[2] * (d_vg - b_vg))
                B_Lvq = np.dot((B_vg[:, X, X] * dddm_vLvg).sum(0), n_qg.T)
                dEdD_sqL[1:] += 2 * pi * w * B_Lvq.transpose((1, 2, 0))
                
            E += w * rgd.integrate(e_g)

        return E, dEdD_sqL

    def calculate_radial(self, rgd, n_Lg, Y_L, dndr_Lg, rnablaY_Lv,
                         m_vLg, dmdr_vLg):
        n_g = np.dot(Y_L, n_Lg)
        m_vg = np.dot(Y_L, m_vLg)
        m_g = (m_vg**2).sum(0)**0.5
        eps = 1e-15
        m_g[m_g < eps] = eps
        n_sg = rgd.empty(2)
        n_sg[:] = n_g
        n_sg[0] += m_g
        n_sg[1] -= m_g
        n_sg *= 0.5

        e_g = rgd.empty()
        dedn_sg = rgd.zeros(2)

        if self.xc.type == 'GGA':
            dmdr_vg = np.dot(Y_L, dmdr_vLg)
        
            a_g = np.dot(Y_L, dndr_Lg)
            b_vg = np.dot(rnablaY_Lv.T, n_Lg)
            
            c_g = (m_vg * dmdr_vg).sum(0) / m_g
            m_vvg = np.dot(rnablaY_Lv.T, m_vLg)
            d_vg = (m_vg * m_vvg).sum(1) / m_g

            sigma_xg = rgd.empty(3)
            sigma_xg[0] = ((b_vg + d_vg)**2).sum(0)
            sigma_xg[1] = ((b_vg + d_vg) * (b_vg - d_vg)).sum(0)
            sigma_xg[2] = ((b_vg - d_vg)**2).sum(0)
            sigma_xg[:, 1:] /= rgd.r_g[1:]**2
            sigma_xg[:, 0] = sigma_xg[:, 1]
            sigma_xg[0] += (a_g + c_g)**2
            sigma_xg[1] += (a_g + c_g) * (a_g - c_g)
            sigma_xg[2] += (a_g - c_g)**2
            sigma_xg *= 0.25

            dedsigma_xg = rgd.zeros(3)

            self.xc.kernel.calculate(e_g, n_sg, dedn_sg, sigma_xg, dedsigma_xg)
        else:
            self.xc.kernel.calculate(e_g, n_sg, dedn_sg)

        dedn_g = 0.5 * dedn_sg.sum(0)
        dedn_sg /= m_g
        dedm_vg = 0.5 * (dedn_sg[0] - dedn_sg[1]) * m_vg
        if self.xc.type == 'GGA':
            dcdm_vg = (dmdr_vg - (m_vg * dmdr_vg).sum(0) * m_vg / m_g**2) / m_g
            dddm_vLvg = rnablaY_Lv.T[:, :, X, X] * m_vg
            dddm_vLvg += m_vvg[:, X] * Y_L[:, X, X]
            dddm_vLvg -= d_vg[:, X, X, :] * m_vg[X, X] * Y_L[:, X, X] / m_g
            dddm_vLvg /= m_g

            return (e_g, dedn_g, dedm_vg,
                    a_g, b_vg, c_g, d_vg, dedsigma_xg, dcdm_vg, dddm_vLvg,
                    m_vg, m_g)

        return e_g, dedn_g, dedm_vg


class NonCollinearLCAOEigensolver(DirectLCAO):
    def calculate_hamiltonian_matrix(self, ham, wfs, kpt, root=-1):

        assert self.has_initialized

        vt_sG = ham.vt_sG
        dH_asp = ham.dH_asp
        H_MM = np.empty((wfs.ksl.mynao, wfs.ksl.nao), wfs.dtype)
        H_sMsM = np.empty((2, wfs.ksl.mynao, 2, wfs.ksl.nao), complex)
        
        wfs.timer.start('Potential matrix')
        self.get_component(wfs, 0, vt_sG, dH_asp, kpt, H_MM)
        H_sMsM[0, :, 0] = H_MM
        H_sMsM[1, :, 1] = H_MM
        self.get_component(wfs, 1, vt_sG, dH_asp, kpt, H_MM)
        H_sMsM[0, :, 1] = H_MM
        H_sMsM[1, :, 0] = H_MM.conj().T
        self.get_component(wfs, 2, vt_sG, dH_asp, kpt, H_MM)
        H_sMsM[0, :, 1] += 1j * H_MM
        H_sMsM[1, :, 0] -= 1j * H_MM.conj().T
        self.get_component(wfs, 3, vt_sG, dH_asp, kpt, H_MM)
        H_sMsM[0, :, 0] += H_MM
        H_sMsM[1, :, 1] -= H_MM
        wfs.timer.stop('Potential matrix')

        ham.gd.comm.sum(H_sMsM)

        H_sMsM[0, :, 0] += wfs.T_qMM[kpt.q]
        H_sMsM[1, :, 1] += wfs.T_qMM[kpt.q]

        wfs.timer.start('Distribute overlap matrix')
        #H_MM = wfs.ksl.distribute_overlap_matrix(H_MM, root)
        wfs.timer.stop('Distribute overlap matrix')
        H_sMsM.shape = (2 * wfs.ksl.mynao, 2 * wfs.ksl.nao)
        return H_sMsM

    def get_component(self, wfs, s, vt_sG, dH_asp, kpt, H_MM):
        wfs.basis_functions.calculate_potential_matrix(vt_sG[s], H_MM, kpt.q)

        # Add atomic contribution
        #
        #           --   a     a  a*
        # H      += >   P    dH  P
        #  mu nu    --   mu i  ij nu j
        #           aij
        #
        wfs.timer.start('Atomic Hamiltonian')
        Mstart = wfs.basis_functions.Mstart
        Mstop = wfs.basis_functions.Mstop
        wfs.timer.stop('Atomic Hamiltonian')
        tri2full(H_MM)
        for a, P_Mi in kpt.P_aMi.items():
            dH_ii = np.asarray(unpack(dH_asp[a][s]), wfs.dtype)
            dHP_iM = np.zeros((dH_ii.shape[1], P_Mi.shape[0]), wfs.dtype)
            # (ATLAS can't handle uninitialized output array)
            gemm(1.0, P_Mi, dH_ii, 0.0, dHP_iM, 'c')
            gemm(1.0, dHP_iM, P_Mi[Mstart:Mstop], 1.0, H_MM)

    def iterate(self, hamiltonian, wfs):
        wfs.timer.start('Noncollinear LCAO eigensolver')

        for kpt in wfs.kpt_u:
            self.iterate_one_k_point(hamiltonian, wfs, kpt)

        wfs.timer.stop('Noncollinear LCAO eigensolver')

    def iterate_one_k_point(self, hamiltonian, wfs, kpt):
        if wfs.bd.comm.size > 1 and wfs.bd.strided:
            raise NotImplementedError

        H_MM = self.calculate_hamiltonian_matrix(hamiltonian, wfs, kpt, root=0)
        S_MM = np.zeros_like(H_MM)
        nao = wfs.ksl.nao
        S_MM.shape = (2, nao, 2, nao)
        for s in range(2):
            S_MM[s, :, s] = wfs.S_qMM[kpt.q]
        S_MM.shape = (2 * nao, 2 * nao)

        if kpt.eps_n is None:
            kpt.eps_n = np.empty(wfs.bd.mynbands)

        diagonalization_string = repr(self.diagonalizer)
        wfs.timer.start(diagonalization_string)
        kpt.C_nsM.shape = (wfs.bd.mynbands, 2 * nao)
        self.diagonalizer.diagonalize(H_MM, kpt.C_nsM, kpt.eps_n, S_MM)
        kpt.C_nsM.shape = (wfs.bd.mynbands, 2, nao)
        wfs.timer.stop(diagonalization_string)


class NonCollinearLCAOWaveFunctions(LCAOWaveFunctions):
    collinear = False
    ncomp = 2

    def set_positions(self, spos_ac):
        LCAOWaveFunctions.set_positions(self, spos_ac)
        for kpt in self.kpt_u:
            kpt.C_nM = None
            kpt.C_nsM = np.empty((self.bd.mynbands, 2, self.ksl.nao), complex)
            
    def add_to_density_from_k_point_with_occupation(self, nt_sG, kpt, f_n):
        rho00_MM = self.ksl.calculate_density_matrix(f_n,
                                                     kpt.C_nsM[:, 0].copy())
        rho01_MM = self.ksl.calculate_density_matrix(
            f_n,
            kpt.C_nsM[:, 0].copy(),
            C2_nM=kpt.C_nsM[:, 1].copy())
        rho11_MM = self.ksl.calculate_density_matrix(f_n,
                                                     kpt.C_nsM[:, 1].copy())
        kpt.rho_sMM = np.empty((4, self.ksl.nao, self.ksl.nao), self.dtype)
        if self.dtype == float:
            kpt.rho_sMM[0] = rho00_MM.real + rho11_MM.real
            kpt.rho_sMM[1] = rho01_MM.real + rho01_MM.T.real
            kpt.rho_sMM[2] = -rho01_MM.imag - rho01_MM.T.imag
            kpt.rho_sMM[3] = rho00_MM.real - rho11_MM.real
        else:
            kpt.rho_sMM[0] = rho00_MM + rho11_MM
            kpt.rho_sMM[1] = rho01_MM + rho01_MM.T.conj()
            kpt.rho_sMM[2] = 1j * (rho01_MM - rho01_MM.T.conj())
            kpt.rho_sMM[3] = rho00_MM - rho11_MM
        for rho_MM, nt_G in zip(kpt.rho_sMM, nt_sG):
            self.basis_functions.construct_density(rho_MM, nt_G, kpt.q)

    def calculate_atomic_density_matrices_k_point(self, D_sii, kpt, a, f_n):
        P_Mi = kpt.P_aMi[a]
        rhoP_Mi = np.zeros_like(P_Mi)
        D0_ii = np.zeros(D_sii[0].shape, self.dtype)
        for rho_MM, D_ii in zip(kpt.rho_sMM, D_sii):
            gemm(1.0, P_Mi, rho_MM, 0.0, rhoP_Mi)
            gemm(1.0, rhoP_Mi, P_Mi.T.conj().copy(), 0.0, D0_ii)
            D_ii += 0.5 * (D0_ii.real + D0_ii.T.real)


class NonCollinearMixer(BaseMixer):
    def mix(self, density):
        BaseMixer.mix(self, density.nt_sG[0],
                      [D_sp[0] for D_sp in density.D_asp.values()])
