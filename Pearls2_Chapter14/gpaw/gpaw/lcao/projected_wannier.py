import numpy as np
import numpy.linalg as la
from ase.units import Hartree

from gpaw.wavefunctions.lcao import add_paw_correction_to_overlap
from gpaw.lcao.overlap import NewTwoCenterIntegrals as TwoCenterIntegrals
from gpaw.lfc import BasisFunctions
from gpaw.utilities import unpack
from gpaw.utilities.tools import dagger, lowdin, tri2full
from gpaw.lcao.tools import basis_subset2


def dots(*args):
    x = args[0]
    for M in args[1:]:
        x = np.dot(x, M)
    return x        


def normalize(U, U2=None, norms=None):
    if norms is None:
        norms2 = np.ones(U.shape[1])
    else:
        norms2 = np.sqrt(norms)

    if U2 is None:
        for n, col in zip(norms2, U.T):
            col *= n / la.norm(col)
    else:
         for n, col1, col2 in zip(norms2, U.T, U2.T):
             norm = np.sqrt(np.vdot(col1, col1) + np.vdot(col2, col2))
             col1 *= n / norm
             col2 *= n / norm
      

def normalize2(C, S):
    C /= np.sqrt(dots(dagger(C), S, C).diagonal())


def get_rot(F_MM, V_oM, L):
    eps_M, U_MM = la.eigh(F_MM)
    indices = eps_M.real.argsort()[-L:] 
    U_Ml = U_MM[:, indices]
    U_Ml /= np.sqrt(dots(U_Ml.T.conj(), F_MM, U_Ml).diagonal())

    U_ow = V_oM.copy()
    U_lw = np.dot(U_Ml.T.conj(), F_MM)
    for col1, col2 in zip(U_ow.T, U_lw.T):
         norm = np.sqrt(np.vdot(col1, col1) + np.vdot(col2, col2))
         col1 /= norm
         col2 /= norm
    return U_ow, U_lw, U_Ml


def condition_number(S):
    eps = la.eigvalsh(S).real
    return eps.max() / eps.min()


def eigvals(H, S):
    return np.sort(la.eigvals(la.solve(S, H)).real)


def get_bfs(calc):
    wfs = calc.wfs
    bfs = BasisFunctions(wfs.gd, [setup.phit_j for setup in wfs.setups],
                         wfs.kd, cut=True)
    bfs.set_positions(calc.atoms.get_scaled_positions() % 1.)
    return bfs


def get_lcao_projections_HSP(calc, bfs=None, spin=0, projectionsonly=True):
    """Some title.

    if projectionsonly is True, return the projections::

      V_qnM = <psi_qn | Phi_qM>

    else, also return the Hamiltonian, overlap, and projector overlaps::

      H_qMM  = <Phi_qM| H |Phi_qM'>
      S_qMM  = <Phi_qM|Phi_qM'>
      P_aqMi = <pt^a_qi|Phi_qM>
    """
    if calc.wfs.kd.comm.size != 1:
        raise NotImplementedError('Parallelization over spin/kpt not '
                                  'implemented yet.')
    spos_ac = calc.atoms.get_scaled_positions() % 1.
    comm = calc.wfs.gd.comm
    nq = len(calc.wfs.kd.ibzk_qc)
    Nk = calc.wfs.kd.nibzkpts
    nao = calc.wfs.setups.nao
    dtype = calc.wfs.dtype
    if bfs is None:
        bfs = get_bfs(calc)
    tci = TwoCenterIntegrals(calc.wfs.gd.cell_cv,
                             calc.wfs.gd.pbc_c,
                             calc.wfs.setups,
                             calc.wfs.kd.ibzk_qc,
                             calc.wfs.kd.gamma)

    # Calculate projector overlaps, and (lower triangle of-) S and T matrices
    S_qMM = np.zeros((nq, nao, nao), dtype)
    T_qMM = np.zeros((nq, nao, nao), dtype)
    P_aqMi = {}

    for a in bfs.my_atom_indices:
        ni = calc.wfs.setups[a].ni
        P_aqMi[a] = np.zeros((nq, nao, ni), dtype)
    tci.calculate(spos_ac, S_qMM, T_qMM, P_aqMi)
    add_paw_correction_to_overlap(calc.wfs.setups, P_aqMi, S_qMM)
    calc.wfs.gd.comm.sum(S_qMM)
    calc.wfs.gd.comm.sum(T_qMM)
    
    
    # Calculate projections
    V_qnM = np.zeros((nq, calc.wfs.bd.nbands, nao), dtype)
    for kpt in calc.wfs.kpt_u:
        if kpt.s != spin:
            continue
        V_nM = V_qnM[kpt.q]
        #bfs.integrate2(kpt.psit_nG[:], V_nM, kpt.q) # all bands to save time
        for n, V_M in enumerate(V_nM): # band-by-band to save memory
            bfs.integrate2(kpt.psit_nG[n][:], V_M, kpt.q)
        for a, P_ni in kpt.P_ani.items():
            dS_ii = calc.wfs.setups[a].dO_ii
            P_Mi = P_aqMi[a][kpt.q]
            V_nM += np.dot(P_ni, np.inner(dS_ii, P_Mi).conj())
    comm.sum(V_qnM)
    if projectionsonly:
        return V_qnM

    # Determine potential matrix
    vt_G = calc.hamiltonian.vt_sG[spin]
    H_qMM = np.zeros((nq, nao, nao), dtype)
    for q, H_MM in enumerate(H_qMM):
        bfs.calculate_potential_matrix(vt_G, H_MM, q)

    # Make Hamiltonian as sum of kinetic (T) and potential (V) matrices
    # and add atomic corrections
    for a, P_qMi in P_aqMi.items():
        dH_ii = unpack(calc.hamiltonian.dH_asp[a][spin])
        for P_Mi, H_MM in zip(P_qMi, H_qMM):
            H_MM += np.dot(P_Mi, np.inner(dH_ii, P_Mi).conj())
    comm.sum(H_qMM)
    H_qMM += T_qMM
    
    # Fill in the upper triangles of H and S
    for H_MM, S_MM in zip(H_qMM, S_qMM):
        tri2full(H_MM)
        tri2full(S_MM)
    H_qMM *= Hartree

    return V_qnM, H_qMM, S_qMM, P_aqMi


def convert_projection_data(symbols, V_qnM, H_qMM, S_qMM, P_aqMi,
                            originaltype='dzp', newtype='sz'):
    Nq = len(H_qMM)
    mask = basis_subset2(symbols, originaltype, newtype)
    mask2 = np.outer(mask, mask)
    NMnew = len(mask.nonzero()[0])
    V_qnM = V_qnM[..., mask]
    H_qMM = H_qMM[:, mask2].reshape(Nq, NMnew, NMnew)
    S_qMM = S_qMM[:, mask2].reshape(Nq, NMnew, NMnew)
    P2_aqMi = {}
    for a, P_qMi in P_aqMi.items():
        P2_aqMi[a] = P_qMi[:, mask, :]
    return V_qnM, H_qMM, S_qMM, P2_aqMi
    

def get_phs(calc, s=0):
    return get_lcao_projections_HSP(calc, bfs=None,
                                    spin=s, projectionsonly=False)


class ProjectedWannierFunctions:
    def __init__(self, projections, h_lcao, s_lcao, eigenvalues, kpoints, 
                 L_k=None, M_k=None, N=None, fixedenergy=None):
        """projections[n,i] = <psi_n|f_i>
           h_lcao[i1, i2] = <f_i1|h|f_i2>
           s_lcao[[i1, i2] = <f_i1|f_i2>
           eps_n: Exact eigenvalues
           L: Number of extra degrees of freedom
           M: Number of states to exactly span
           N: Total number of bands in the calculation
           
           Methods:
           -- get_hamiltonian_and_overlap_matrix --
           will return the hamiltonian and identity operator
           in the projected wannier function basis. 
           The following steps are performed:
            
           1) calculate_edf       -> self.b_il
           2) calculate_rotations -> self.Uo_mi and self.Uu_li
           3) calculate_overlaps  -> self.S_ii
           4) calculate_hamiltonian_matrix -> self.H_ii

           -- get_eigenvalues --
           gives the eigenvalues of of the hamiltonian in the
           projected wannier function basis.

           -- indices --
           i localized function index
           n eigenstate index
           l edf index
           m fixed eigenstate index
           k k-point index
           """
         
        self.eps_kn = eigenvalues
        self.ibzk_kc = kpoints
        self.nk = len(self.ibzk_kc)
        self.V_kni = projections   #<psi_n1|f_i1>
        self.dtype = self.V_kni.dtype
        self.Nw = self.V_kni.shape[2]
        self.s_lcao_kii = s_lcao #F_ii[i1,i2] = <f_i1|f_i2>
        self.h_lcao_kii = h_lcao

        if N is None:
            N = self.V_kni.shape[1]
        self.N = N
        
        if fixedenergy is None:
            raise NotImplementedError,'Only fixedenergy is implemented for now'
        else:
            self.fixedenergy = fixedenergy
            self.M_k = [sum(eps_n <= fixedenergy) for eps_n in self.eps_kn]
            self.L_k = [self.Nw - M for M in self.M_k]
            print("fixedenergy =", self.fixedenergy)

        print('N =', self.N)
        print('skpt_kc = ')
        print(self.ibzk_kc)
        print('M_k =', self.M_k)
        print('L_k =', self.L_k)
        print('Nw =', self.Nw)

    def get_hamiltonian_and_overlap_matrix(self, useibl=True):
        self.calculate_edf(useibl=useibl)
        self.calculate_rotations()
        self.calculate_overlaps()
        self.calculate_hamiltonian_matrix(useibl=useibl)
        return self.H_kii, self.S_kii

    def calculate_edf(self, useibl=True):
        """Calculate the coefficients b_il in the expansion of the EDF.

        ``|phi_l> = sum_i b_il |f^u_i>``, in terms of ``|f^u_i> = P^u|f_i>``.

        To use the infinite band limit set useibl=True.
        N is the total number of bands to use.
        """
        
        for k, L in enumerate(self.L_k):
            if L==0:
                assert L!=0, 'L_k=0 for k=%i. Not implemented' % k
        
        self.Vo_kni = [V_ni[:M] for V_ni, M in zip(self.V_kni, self.M_k)]
        
        self.Fo_kii = np.asarray([np.dot(dagger(Vo_ni), Vo_ni) 
                                  for Vo_ni in self.Vo_kni])
        
        if useibl:
            self.Fu_kii = self.s_lcao_kii - self.Fo_kii
        else:
            self.Vu_kni = [V_ni[M:self.N] 
                           for V_ni, M in zip(self.V_kni, self.M_k)]
            self.Fu_kii = np.asarray([np.dot(dagger(Vu_ni), Vu_ni) 
                                     for Vu_ni in self.Vu_kni])
        self.b_kil = [] 
        for Fu_ii, L in zip(self.Fu_kii, self.L_k):
            b_i, b_ii = la.eigh(Fu_ii)
            ls = b_i.real.argsort()[-L:] 
            b_il = b_ii[:, ls] #pick out the eigenvec with largest eigenvals.
            normalize2(b_il, Fu_ii) #normalize the EDF: <phi_l|phi_l> = 1
            self.b_kil.append(b_il)

    def calculate_rotations(self):
        Uo_kni = [Vo_ni.copy() for Vo_ni in self.Vo_kni]
        Uu_kli = [np.dot(dagger(b_il), Fu_ii) 
                  for b_il, Fu_ii in zip(self.b_kil, self.Fu_kii)]
        #Normalize such that <omega_i|omega_i> = <f_i|f_i>
        for Uo_ni, Uu_li, s_ii in zip(Uo_kni, Uu_kli, self.s_lcao_kii):
            normalize(Uo_ni, Uu_li, s_ii.diagonal())
        self.Uo_kni = Uo_kni
        self.Uu_kli = Uu_kli

    def calculate_overlaps(self):
        Wo_kii = [np.dot(dagger(Uo_ni), Uo_ni) for Uo_ni in self.Uo_kni]
        Wu_kii = [np.dot(dagger(Uu_li), Uu_li) for Uu_li in self.Uu_kli]
        #Wu_kii = [dots(dagger(Uu_li), dagger(b_il), Fu_ii, b_il, Uu_li) 
        #for Uu_li, b_il, Fu_ii in zip(self.Uu_kli, self.b_kil, self.Fu_kii)]
        self.S_kii = np.asarray(Wo_kii) + np.asarray(Wu_kii)

    def get_condition_number(self):
        eigs_kn = [la.eigvalsh(S_ii) for S_ii in self.S_kii]
        return np.asarray([condition_number(S) for S in self.S_kii])

    def calculate_hamiltonian_matrix(self, useibl=True):
        """Calculate H_kij = H^o_i(k)j(k) + H^u_i(k)j(k)
           i(k): Bloch sum of omega_i
        """

        epso_kn = [eps_n[:M] for eps_n, M in zip(self.eps_kn, self.M_k)]
        self.Ho_kii = np.asarray([np.dot(dagger(Uo_ni) * epso_n, Uo_ni) 
                                  for Uo_ni, epso_n in zip(self.Uo_kni, 
                                                           epso_kn)])

        if self.h_lcao_kii is not None and useibl:
            print("Using h_lcao and infinite band limit")
            Vo_kni = self.Vo_kni
            Huf_kii = [h_lcao_ii - np.dot(dagger(Vo_ni) * epso_n, Vo_ni)
                       for h_lcao_ii, Vo_ni, epso_n in zip(self.h_lcao_kii, 
                                                           self.Vo_kni, 
                                                           epso_kn)]
            self.Huf_kii = np.asarray(Huf_kii)
        else:
            print("Using finite band limit (not using h_lcao)")
            epsu_kn = [eps_n[M:self.N] 
                       for eps_n, M in zip(self.eps_kn, self.M_k)]
            Huf_kii = [np.dot(dagger(Vu_ni) * epsu_n, Vu_ni) 
                       for Vu_ni, epsu_n in zip(self.Vu_kni, epsu_kn)]
            self.Huf_kii = np.asarray(Huf_kii)

        Hu_kii = [dots(dagger(Uu_li), dagger(b_il), Huf_ii, b_il, Uu_li)
                  for Uu_li, b_il, Huf_ii in zip(self.Uu_kli, self.b_kil,
                                                 self.Huf_kii)]
        self.Hu_kii = np.asarray(Hu_kii)
        self.H_kii = self.Ho_kii + self.Hu_kii

    def get_eigenvalues(self):
        return np.asarray([eigvals(H, S)
                           for H, S in zip(self.H_kii, self.S_kii)])

    def get_lcao_eigenvalues(self):
        return np.asarray([eigvals(H, S)
                           for H, S in zip(self.h_lcao_kii, self.s_lcao_kii)])

    def get_norm_of_projection(self):
        norm_kn = np.zeros((self.nk, self.N))
        Sinv_kii = np.asarray([la.inv(S_ii) for S_ii in self.S_kii])

        normo_kn = np.asarray([dots(Uo_ni, Sinv_ii, dagger(Uo_ni)).diagonal()
                    for Uo_ni, Sinv_ii in zip(self.Uo_kni, Sinv_kii)])
        
        Vu_kni = np.asarray([V_ni[M:self.N] 
                             for V_ni, M in zip(self.V_kni, self.M_k)])

        Pu_kni = [dots(Vu_ni, b_il, Uu_li) 
                for Vu_ni, b_il, Uu_li in zip(Vu_kni, self.b_kil, self.Uu_kli)]

        normu_kn = np.asarray([dots(Pu_ni, Sinv_ii, dagger(Pu_ni)).diagonal()
                    for Pu_ni, Sinv_ii in zip(Pu_kni, Sinv_kii)])

        return np.hstack((normo_kn, normu_kn))
        

    def calculate_functions(self, calc, basis, k=0):
        from gpaw.io import FileReference
        psit_nG = calc.wfs.kpt_u[k].psit_nG
        atoms = calc.get_atoms()
        Uo_ni = self.Uo_kni[k]
        tarinstance = isinstance(psit_nG, FileReference)
        if tarinstance:
            psit_nG = np.asarray([psit_nG[i] for i in range(self.M_k[k])])

        basis_functions = get_bfs(calc)
        b_il, Uu_li, Vo_ni = self.b_kil[k], self.Uu_kli[k], self.Vo_kni[k]
        a_iG = -np.tensordot(Vo_ni, psit_nG, axes=[0, 0])# a_iG = -fo_iG
        #self.fo_iG = -a_iG.copy()#f projected onto the occupied subspace
        C_iM = np.identity(self.Nw, dtype=self.dtype)
        basis_functions.lcao_to_grid(C_iM, a_iG, q=-1) # a_iG=fu_iG=f_iG-fo_iG
        #self.f_iG = self.fo_iG + a_iG # check
        a_iG = np.tensordot(b_il, a_iG, axes=[0, 0])   # a_iG = EDF
        a_iG = np.tensordot(Uu_li, a_iG, axes=[0, 0])  # a_iG = wu_iG
        a_iG+=np.tensordot(Uo_ni, psit_nG, axes=[0, 0])# ai_G = wu_iG+wo_iG
        self.w_iG = a_iG


    def get_mlwf_initial_guess(self):
        """calculate initial guess for maximally localized 
        wannier functions. Does not work for the infinite band limit.
        cu_nl: rotation coefficents of unoccupied states
        U_ii: rotation matrix of eigenstates and edf.
        """
        Vu_ni = self.Vu_ni[self.M: self.N]
        cu_nl = np.dot(Vu_ni, self.b_il)
        U_ii = np.vstack((self.Uo_ni, self.Uu_li))
        lowdin(U_ii)
        return U_ii, cu_nl
