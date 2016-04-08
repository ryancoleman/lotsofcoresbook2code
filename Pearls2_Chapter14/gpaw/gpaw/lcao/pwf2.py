import numpy as np

from ase.units import Hartree
from gpaw.aseinterface import GPAW
from gpaw.lcao.overlap import NewTwoCenterIntegrals
from gpaw.utilities import unpack
from gpaw.utilities.tools import tri2full, lowdin
from gpaw.lcao.tools import basis_subset2, get_bfi2
from gpaw.coulomb import get_vxc as get_ks_xc
from gpaw.utilities.blas import r2k, gemm
from gpaw.lcao.projected_wannier import dots, condition_number, eigvals, \
     get_bfs, get_lcao_projections_HSP


def get_rot(F_MM, V_oM, L):
    eps_M, U_MM = np.linalg.eigh(F_MM)
    indices = eps_M.real.argsort()[-L:] 
    U_Ml = U_MM[:, indices]
    U_Ml /= np.sqrt(dots(U_Ml.T.conj(), F_MM, U_Ml).diagonal())

    U_ow = V_oM.copy()
    U_lw = np.dot(U_Ml.T.conj(), F_MM)
    for col1, col2 in zip(U_ow.T, U_lw.T):
        norm = np.linalg.norm(np.hstack((col1, col2)))
        col1 /= norm
        col2 /= norm
    return U_ow, U_lw, U_Ml
    

def get_lcao_xc(calc, P_aqMi, bfs=None, spin=0):
    nq = len(calc.wfs.kd.ibzk_qc)
    nao = calc.wfs.setups.nao
    dtype = calc.wfs.dtype
    if bfs is None:
        bfs = get_bfs(calc)
    
    if calc.density.nt_sg is None:
        calc.density.interpolate_pseudo_density()
    nt_sg = calc.density.nt_sg
    vxct_sg = calc.density.finegd.zeros(calc.wfs.nspins)
    calc.hamiltonian.xc.calculate(calc.density.finegd, nt_sg, vxct_sg)
    vxct_G = calc.wfs.gd.zeros()
    calc.hamiltonian.restrict(vxct_sg[spin], vxct_G)
    Vxc_qMM = np.zeros((nq, nao, nao), dtype)
    for q, Vxc_MM in enumerate(Vxc_qMM):
        bfs.calculate_potential_matrix(vxct_G, Vxc_MM, q)
        tri2full(Vxc_MM, 'L')

    # Add atomic PAW corrections
    for a, P_qMi in P_aqMi.items():
        D_sp = calc.density.D_asp[a][:]
        H_sp = np.zeros_like(D_sp)
        calc.hamiltonian.xc.calculate_paw_correction(calc.wfs.setups[a],
                                                     D_sp, H_sp) 
        H_ii = unpack(H_sp[spin])
        for Vxc_MM, P_Mi in zip(Vxc_qMM, P_qMi):
            Vxc_MM += dots(P_Mi, H_ii, P_Mi.T.conj())
    return Vxc_qMM * Hartree


def get_xc2(calc, w_wG, P_awi, spin=0):
    if calc.density.nt_sg is None:
        calc.density.interpolate_pseudo_density()
    nt_g = calc.density.nt_sg[spin]
    vxct_g = calc.density.finegd.zeros()
    calc.hamiltonian.xc.get_energy_and_potential(nt_g, vxct_g)
    vxct_G = calc.wfs.gd.empty()
    calc.hamiltonian.restrict(vxct_g, vxct_G)

    # Integrate pseudo part
    Nw = len(w_wG)
    xc_ww = np.empty((Nw, Nw))
    r2k(.5 * calc.wfs.gd.dv, w_wG, vxct_G * w_wG, .0, xc_ww)
    tri2full(xc_ww, 'L')
    
    # Add atomic PAW corrections
    for a, P_wi in P_awi.items():
        D_sp = calc.density.D_asp[a][:]
        H_sp = np.zeros_like(D_sp)
        calc.wfs.setups[a].xc_correction.calculate_energy_and_derivatives(
            D_sp, H_sp)
        H_ii = unpack(H_sp[spin])
        xc_ww += dots(P_wi, H_ii, P_wi.T.conj())
    return xc_ww * Hartree


class ProjectedWannierFunctionsFBL:
    """PWF in the finite band limit.

    ::
    
                --N              
        |w_w> = >    |psi_n> U_nw
                --n=1            
    """
    def __init__(self, V_nM, No, ortho=False):
        Nw = V_nM.shape[1]
        assert No <= Nw
        
        V_oM, V_uM = V_nM[:No], V_nM[No:]
        F_MM = np.dot(V_uM.T.conj(), V_uM)
        U_ow, U_lw, U_Ml = get_rot(F_MM, V_oM, Nw - No)
        self.U_nw = np.vstack((U_ow, dots(V_uM, U_Ml, U_lw)))

        # stop here ?? XXX
        self.S_ww = self.rotate_matrix(np.ones(1))
        if ortho:
            lowdin(self.U_nw, self.S_ww)
            self.S_ww = np.identity(Nw)
        self.norms_n = np.dot(self.U_nw, np.linalg.solve(
            self.S_ww, self.U_nw.T.conj())).diagonal()

    def rotate_matrix(self, A_nn):
        if A_nn.ndim == 1:
            return np.dot(self.U_nw.T.conj() * A_nn, self.U_nw)
        else:
            return dots(self.U_nw.T.conj(), A_nn, self.U_nw)

    def rotate_projections(self, P_ani):
        P_awi = {}
        for a, P_ni in P_ani.items():
            P_awi[a] = np.tensordot(self.U_nw, P_ni, axes=[[0], [0]])
        return P_awi

    def rotate_function(self, psit_nG):
        return np.tensordot(self.U_nw, psit_nG, axes=[[0], [0]])


class ProjectedWannierFunctionsIBL:
    """PWF in the infinite band limit.

    ::
    
                --No               --Nw
        |w_w> = >   |psi_o> U_ow + >   |f_M> U_Mw
                --o=1              --M=1
    """
    def __init__(self, V_nM, S_MM, No, lcaoindices=None):
        Nw = V_nM.shape[1]
        assert No <= Nw
        self.V_oM, V_uM = V_nM[:No], V_nM[No:]
        
        F_MM = S_MM - np.dot(self.V_oM.T.conj(), self.V_oM)
        U_ow, U_lw, U_Ml = get_rot(F_MM, self.V_oM, Nw - No)
        self.U_Mw = np.dot(U_Ml, U_lw)
        self.U_ow = U_ow - np.dot(self.V_oM, self.U_Mw)
        if lcaoindices is not None:
            for i in lcaoindices:
                self.U_ow[:, i] = 0.0
                self.U_Mw[:, i] = 0.0
                self.U_Mw[i, i] = 1.0

        # stop here ?? XXX
        self.S_ww = self.rotate_matrix(np.ones(1), S_MM)
        P_uw = np.dot(V_uM, self.U_Mw)
        self.norms_n = np.hstack((
           np.dot(U_ow, np.linalg.solve(self.S_ww, U_ow.T.conj())).diagonal(),
           np.dot(P_uw, np.linalg.solve(self.S_ww, P_uw.T.conj())).diagonal()))

    def rotate_matrix(self, A_o, A_MM):
        assert A_o.ndim == 1
        A_ww = dots(self.U_ow.T.conj() * A_o, self.V_oM, self.U_Mw)
        A_ww += np.conj(A_ww.T)
        A_ww += np.dot(self.U_ow.T.conj() * A_o, self.U_ow)
        A_ww += dots(self.U_Mw.T.conj(), A_MM, self.U_Mw)
        return A_ww

    def rotate_projections(self, P_aoi, P_aMi, indices=None):
        if indices is None:
            U_ow = self.U_ow
            U_Mw = self.U_Mw
        else:
            U_ow = self.U_ow[:, indices]
            U_Mw = self.U_Mw[:, indices]
        P_awi = {}
        for a, P_oi in P_aoi.items():
            P_awi[a] = np.tensordot(U_Mw, P_aMi[a], axes=[[0], [0]])
            if len(U_ow) > 0:
                P_awi[a] += np.tensordot(U_ow, P_oi, axes=[[0], [0]])
        return P_awi

    def rotate_function(self, psit_oG, bfs, q=-1, indices=None):
        if indices is None:
            U_ow = self.U_ow
            U_Mw = self.U_Mw
        else:
            U_ow = self.U_ow[:, indices]
            U_Mw = self.U_Mw[:, indices]
        w_wG = np.zeros((U_ow.shape[1],) + psit_oG.shape[1:])
        if len(U_ow) > 0:
            gemm(1., psit_oG, U_ow.T.copy(), 0., w_wG)
        bfs.lcao_to_grid(U_Mw.T.copy(), w_wG, q)
        return w_wG


class PWFplusLCAO(ProjectedWannierFunctionsIBL):
    def __init__(self, V_nM, S_MM, No, pwfmask, lcaoindices=None):
        Nw = V_nM.shape[1]
        self.V_oM = V_nM[:No]
        dtype = V_nM.dtype
        
        # Do PWF optimization for pwfbasis submatrix only!
        Npwf = len(pwfmask.nonzero()[0])
        pwfmask2 = np.outer(pwfmask, pwfmask)
        s_MM = S_MM[pwfmask2].reshape(Npwf, Npwf)
        v_oM = self.V_oM[:, pwfmask]
        f_MM = s_MM - np.dot(v_oM.T.conj(), v_oM)
        nw = len(s_MM)
        assert No <= nw
        
        u_ow, u_lw, u_Ml = get_rot(f_MM, v_oM, nw - No)
        u_Mw = np.dot(u_Ml, u_lw)
        u_ow = u_ow - np.dot(v_oM, u_Mw)

        # Determine U for full lcao basis
        self.U_ow = np.zeros((No, Nw), dtype)
        for U_w, u_w in zip(self.U_ow, u_ow):
            np.place(U_w, pwfmask, u_w)
        self.U_Mw = np.identity(Nw, dtype)
        np.place(self.U_Mw, pwfmask2, u_Mw.flat)

        if lcaoindices is not None:
            for i in lcaoindices:
                self.U_ow[:, i] = 0.0
                self.U_Mw[:, i] = 0.0
                self.U_Mw[i, i] = 1.0

        self.S_ww = self.rotate_matrix(np.ones(1), S_MM)
        self.norms_n = None

def set_lcaoatoms(calc, pwf, lcaoatoms):
    ind = get_bfi(calc, lcaoatoms)
    for i in ind:
        pwf.U_ow[:, i] = 0.0
        pwf.U_Mw[:, i] = 0.0
        pwf_U_Mw[i, i] = 1.0

class PWF2:
    def __init__(self, gpwfilename, fixedenergy=0., spin=0, ibl=True,
                 basis='sz', zero_fermi=False, pwfbasis=None, lcaoatoms=None,
                 projection_data=None):
        calc = GPAW(gpwfilename, txt=None, basis=basis)
        assert calc.wfs.gd.comm.size == 1
        assert calc.wfs.kd.comm.size == 1
        assert calc.wfs.band_comm.size == 1
        if zero_fermi:
            try:
                Ef = calc.get_fermi_level()
            except NotImplementedError:
                Ef = calc.get_homo_lumo().mean()
        else:
            Ef = 0.0

        self.ibzk_kc = calc.get_ibz_k_points()
        self.nk = len(self.ibzk_kc)
        self.eps_kn = [calc.get_eigenvalues(kpt=q, spin=spin) - Ef
                       for q in range(self.nk)]
        self.M_k = [sum(eps_n <= fixedenergy) for eps_n in self.eps_kn]
        print('Fixed states:', self.M_k) 
        self.calc = calc
        self.dtype = self.calc.wfs.dtype
        self.spin = spin
        self.ibl = ibl
        self.pwf_q = []
        self.norms_qn = []
        self.S_qww = []
        self.H_qww = []

        if ibl:
            if pwfbasis is not None:
                pwfmask = basis_subset2(calc.atoms.get_chemical_symbols(),
                                       basis, pwfbasis)
            if lcaoatoms is not None:
                lcaoindices = get_bfi2(calc.atoms.get_chemical_symbols(),
                                       basis,
                                       lcaoatoms)
            else:
                lcaoindices = None
            self.bfs = get_bfs(calc)
            if projection_data is None:
                V_qnM, H_qMM, S_qMM, self.P_aqMi = get_lcao_projections_HSP(
                    calc, bfs=self.bfs, spin=spin, projectionsonly=False)
            else:
                V_qnM, H_qMM, S_qMM, self.P_aqMi = projection_data
            H_qMM -= Ef * S_qMM
            for q, M in enumerate(self.M_k):
                if pwfbasis is None:
                    pwf = ProjectedWannierFunctionsIBL(V_qnM[q], S_qMM[q], M,
                                                       lcaoindices)
                else:
                    pwf = PWFplusLCAO(V_qnM[q], S_qMM[q], M, pwfmask,
                                      lcaoindices)
                self.pwf_q.append(pwf)
                self.norms_qn.append(pwf.norms_n)
                self.S_qww.append(pwf.S_ww)
                self.H_qww.append(pwf.rotate_matrix(self.eps_kn[q][:M],
                                                    H_qMM[q]))
        else:
            if projection_data is None:
                V_qnM = get_lcao_projections_HSP(calc, spin=spin)
            else:
                V_qnM = projection_data
            for q, M in enumerate(self.M_k):
                pwf = ProjectedWannierFunctionsFBL(V_qnM[q], M, ortho=False)
                self.pwf_q.append(pwf)
                self.norms_qn.append(pwf.norms_n)
                self.S_qww.append(pwf.S_ww)
                self.H_qww.append(pwf.rotate_matrix(self.eps_kn[q]))

        for S in self.S_qww:
            print('Condition number: %0.1e' % condition_number(S))

    def get_hamiltonian(self, q=0, indices=None):
        if indices is None:
            return self.H_qww[q]
        else:
            return self.H_qww[q].take(indices, 0).take(indices, 1)

    def get_overlap(self, q=0, indices=None):
        if indices is None:
            return self.S_qww[q]
        else:
            return self.S_qww[q].take(indices, 0).take(indices, 1)

    def get_projections(self, q=0, indices=None):
        kpt = self.calc.wfs.kpt_u[self.spin * self.nk + q]
        if not hasattr(self, 'P_awi'):
            if self.ibl:
                M = self.M_k[q]
                self.P_awi = self.pwf_q[q].rotate_projections(
                    dict([(a, P_ni[:M]) for a, P_ni in kpt.P_ani.items()]),
                    dict([(a, P_qMi[q]) for a, P_qMi in self.P_aqMi.items()]),
                    indices)
            else:
                self.P_awi = pwf.rotate_projections(kpt.P_ani, indices)
        return self.P_awi

    def get_orbitals(self, q=0, indices=None):
        self.calc.wfs.initialize_wave_functions_from_restart_file()
        kpt = self.calc.wfs.kpt_u[self.spin * self.nk + q]
        if not hasattr(self, 'w_wG'):
            if self.ibl:
                self.w_wG = self.pwf_q[q].rotate_function(
                    kpt.psit_nG[:self.M_k[q]], self.bfs, q, indices)
            else:
                self.w_wG = self.pwf_q[q].rotate_function(
                    kpt.psit_nG, indices)
        return self.w_wG

    def get_Fcore(self, q=0, indices=None):
        if indices is None:
            Fcore_ww = np.zeros_like(self.H_qww[q])
        else:
            Fcore_ww = np.zeros((len(indices), len(indices)))
        for a, P_wi in self.get_projections(q, indices).items():
            X_ii = unpack(self.calc.wfs.setups[a].X_p)
            Fcore_ww -= dots(P_wi, X_ii, P_wi.T.conj())
        return Fcore_ww * Hartree

    def get_eigs(self, q=0):
        return eigvals(self.H_qww[q], self.S_ww[q])

    def get_condition_number(self, q=0):
        return condition_number(self.S_qww[q])

    def get_xc(self, q=0, indices=None):
        #self.calc.density.ghat.set_positions(
        #    self.calc.atoms.get_scaled_positions() % 1.)
        #self.calc.hamiltonian.poisson.initialize()
        if self.ibl:
            return get_xc2(self.calc, self.get_orbitals(q, indices),
                           self.get_projections(q, indices), self.spin)
        else:
            return self.pwf_q[q].rotate_matrix(get_ks_xc(self.calc,
                                                         spin=self.spin))


class LCAOwrap:
    def __init__(self, calc, spin=0):
        assert calc.wfs.gd.comm.size == 1
        assert calc.wfs.kd.comm.size == 1
        assert calc.wfs.band_comm.size == 1
        
        from gpaw.lcao.tools import get_lcao_hamiltonian
        H_skMM, S_kMM = get_lcao_hamiltonian(calc)

        self.calc = calc
        self.dtype = calc.wfs.dtype
        self.spin = spin
        self.H_qww = H_skMM[spin]
        self.S_qww = S_kMM
        self.P_aqwi = calc.wfs.P_aqMi
        self.Nw = self.S_qww.shape[-1]
        
        for S in self.S_qww:
            print('Condition number: %0.1e' % condition_number(S))
    
    def get_hamiltonian(self, q=0, indices=None):
        if indices is None:
            return self.H_qww[q]
        else:
            return self.H_qww[q].take(indices, 0).take(indices, 1)

    def get_overlap(self, q=0, indices=None):
        if indices is None:
            return self.S_qww[q]
        else:
            return self.S_qww[q].take(indices, 0).take(indices, 1)

    def get_projections(self, q=0, indices=None):
        if indices is None:
            return dict([(a, P_qwi[q]) for a, P_qwi in self.P_aqwi.items()])
        else:
            return dict([(a, P_qwi[q].take(indices, 0))
                         for a, P_qwi in self.P_aqwi.items()])

    def get_orbitals(self, q=-1, indices=None):
        assert q == -1
        if indices is None:
            indices = range(self.Nw)
        Ni = len(indices)
        C_wM = np.zeros((Ni, self.Nw), self.dtype)
        for i, C_M in zip(indices, C_wM):
            C_M[i] = 1.0
        w_wG = self.calc.wfs.gd.zeros(Ni, dtype=self.dtype)
        self.calc.wfs.basis_functions.lcao_to_grid(C_wM, w_wG, q=-1)
        return w_wG

    def get_Fcore(self, q=0, indices=None):
        if indices is None:
            Fcore_ww = np.zeros_like(self.H_qww[q])
        else:
            Fcore_ww = np.zeros((len(indices), len(indices)))
        for a, P_wi in self.get_projections(q, indices).items():
            if self.calc.wfs.setups[a].type != 'ghost':
                X_ii = unpack(self.calc.wfs.setups[a].X_p)
                Fcore_ww -= dots(P_wi, X_ii, P_wi.T.conj())
        return Fcore_ww * Hartree

    def get_xc(self, q=0, indices=None):
        if not hasattr(self, 'Vxc_qww'):
            self.Vxc_qww = get_lcao_xc(self.calc, self.P_aqwi,
                                       bfs=self.calc.wfs.basis_functions,
                                       spin=self.spin)
        if indices is None:
            return self.Vxc_qww[q]
        else:
            return self.Vxc_qww[q].take(indices, 0).take(indices, 1)
