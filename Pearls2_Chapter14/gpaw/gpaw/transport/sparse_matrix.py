import numpy as np

from ase.utils.timing import Timer

from gpaw.utilities import unpack
from gpaw.mpi import world, rank
from gpaw.utilities.blas import gemm
from gpaw.utilities.lapack import inverse_general
from gpaw.transport.tools import get_matrix_index, collect_lead_mat, dot
import copy
import _gpaw

class Banded_Sparse_HSD:
    #for lead's hamiltonian, overlap, and density matrix
    def __init__(self, dtype, ns, npk, index=None):
        self.band_index = index
        self.dtype = dtype
        self.H = []
        self.S = []
        self.D = []
        self.ns = ns
        self.npk = npk
        self.s = 0
        self.pk = 0
        for s in range(ns):
            self.H.append([])
            self.D.append([])
            for k in range(npk):
                self.H[s].append([])
                self.D[s].append([])
        for k in range(npk):
            self.S.append([])        

    def reset(self, s, pk, mat, flag='S', init=False):
        assert mat.dtype == self.dtype
        if flag == 'S':
            spar = self.S
        elif flag == 'H':
            spar = self.H[s]
        elif flag == 'D':
            spar = self.D[s]
        if not init:
            spar[pk].reset(mat)            
        elif self.band_index is not None:
            spar[pk] = Banded_Sparse_Matrix(self.dtype, mat, self.band_index)
        else:
            spar[pk] = Banded_Sparse_Matrix(self.dtype, mat)
            self.band_index = spar[pk].band_index
       
class Banded_Sparse_Matrix:
    def __init__(self, dtype, mat=None, band_index=None, tol=1e-9):
        self.tol = tol
        self.dtype = dtype
        self.band_index = band_index
        if mat is not None:
            if band_index is None:
                self.initialize(mat)
            else:
                self.reset(mat)
           
    def initialize(self, mat):
        # the indexing way needs mat[0][-1] = 0,otherwise will recover a
        # unsymmetric full matrix
        assert self.dtype == mat.dtype
        #assert mat[0][-1] < self.tol
        dim = mat.shape[-1]
        ku = -1
        kl = -1
        mat_sum = np.sum(abs(mat))
        spar_sum = 0
        while abs(mat_sum - spar_sum) > self.tol * 10:
            #ku += 1
            #kl += 1
            #ud_sum = 1
            #dd_sum = 1
            #while(ud_sum > self.tol):
            #    ku += 1
            #    ud_sum = np.sum(np.diag(abs(mat), ku))
            #while(dd_sum > self.tol):
            #    kl += 1
            #    dd_sum = np.sum(np.diag(abs(mat), -kl))
            #ku -= 1
            #kl -= 1
            ku = dim
            kl = dim
   
            # storage in the tranpose, bacause column major order for zgbsv_ function
            length = (kl + ku + 1) * dim - kl * (kl + 1) / 2. - \
                                                ku * (ku + 1) / 2.
            
            self.spar = np.zeros([length], self.dtype)
                
            index1 = []
            index2 = []
            index0 = np.zeros((dim, 2 * kl + ku + 1), int)
           
            n = 0
            for i in range(kl, -1, -1):
                for j in range(dim - i):
                    index1.append(i + j)
                    index2.append(j)
                    index0[i + j, 2 * kl - i] = n
                    n += 1
        
            for i in range(1, ku + 1):
                for j in range(dim - i):
                    index1.append(j)
                    index2.append(j + i)
                    index0[j, 2 * kl + i] = n
                    n += 1
            
            index1 = np.array(index1)        
            index2 = np.array(index2)
            
            self.band_index = (kl, ku, index0, index1, index2)
            self.spar = mat[index1, index2]
            spar_sum = np.sum(abs(self.recover()))

    def test1(self, n1, n2):
        index1 ,index2 = self.band_index[-2:]
        for i in range(len(index1)):
            if index1[i] == n1 and index2[i] == n2:
                print(i)
                
    def recover(self):
        index0, index1, index2 = self.band_index[-3:]
        dim = index0.shape[0]
        mat = np.zeros([dim, dim], self.dtype)
        mat[index1, index2] = self.spar
        return mat
 
    def reset(self, mat):
        index1, index2 = self.band_index[-2:]
        assert self.dtype == mat.dtype
        self.spar = mat[index1, index2]

    def reset_from_others(self, bds_mm1, bds_mm2, c1, c2):
        assert self.dtype == complex
        self.spar = c1 * bds_mm1.spar + c2 * bds_mm2.spar 
            
    def reset_minus(self, mat, full=False):
        assert self.dtype == complex
        index1, index2 = self.band_index[-2:]
        if full:
            self.spar -= mat[index1, index2]
        else:
            self.spar -= mat.recover()[index1, index2]
    
    def reset_plus(self, mat, full=False):
        assert self.dtype == complex
        index1, index2 = self.band_index[-2:]
        if full:
            self.spar += mat[index1, index2]
        else:
            self.spar += mat.recover()[index1, index2]           

    def test_inv_speed(self):
        full_mat = self.recover()
        timer = Timer()
        timer.start('full_numpy')
        tmp0 = np.linalg.inv(full_mat)
        timer.stop('full_numpy')
        
        timer.start('full_lapack')
        inverse_general(full_mat)
        timer.stop('full_lapack')
        
        timer.start('sparse_lapack')
        self.inv()
        timer.stop('sparse_lapack')
        
        times = []
        methods = ['full_numpy', 'full_lapack', 'sparse_lapack']
        for name in methods:
            time = timer.timers[name,]
            print(name, time)
            times.append(time)
        
        mintime = np.min(times)
        self.inv_method = methods[np.argmin(times)]
        print('mintime', mintime)
                
    def inv(self):
        #kl, ku, index0 = self.band_index[:3]
        #dim = index0.shape[0]
        #inv_mat = np.eye(dim, dtype=complex)
        #ldab = 2*kl + ku + 1
        #source_mat = self.spar[index0]
        #assert source_mat.flags.contiguous
        #info = _gpaw.linear_solve_band(source_mat, inv_mat, kl, ku)            
        #return inv_mat
        return np.linalg.inv(self.recover()).copy()
       
class Tp_Sparse_HSD:
    def __init__(self, dtype, ns, npk, ll_index, ex=True):
        self.dtype = dtype
        self.ll_index = ll_index
        self.extended = ex
        self.H = []
        self.S = []
        self.D = []
        self.G = []
        self.ns = ns
        self.npk = npk
        self.s = 0
        self.pk = 0
        self.band_indices = None
        for s in range(ns):
            self.H.append([])
            self.D.append([])
            for k in range(npk):
                self.H[s].append([])
                self.D[s].append([])
        for k in range(npk):
            self.S.append([])
        self.G = Tp_Sparse_Matrix(complex, self.ll_index,
                                                    None, None, self.extended)
    
    def reset(self, s, pk, mat, flag='S', init=False):
        if flag == 'S':
            spar = self.S
        elif flag == 'H':
            spar = self.H[s]
        elif flag == 'D':
            spar = self.D[s]
        if not init:
            spar[pk].reset(mat)
        elif self.band_indices is None:
            spar[pk] = Tp_Sparse_Matrix(self.dtype, self.ll_index, mat,
                                                          None, self.extended)
            self.band_indices = spar[pk].band_indices
        else:
            spar[pk] = Tp_Sparse_Matrix(self.dtype, self.ll_index, mat,
                                             self.band_indices, self.extended)

    def append_lead_as_buffer(self, lead_hsd, lead_couple_hsd, ex_index, tp=None):
        assert self.extended == True
        clm = collect_lead_mat
        if tp is not None:
            tp.log('append_lead_as_buffer(), npk : {0}  ns : {1}'.format(self.npk, self.ns))
        for pk in range(self.npk):
            if tp is not None:
                tp.log('append_lead_as_buffer(), pk : {0}'.format(pk))
            diag_h, upc_h, dwnc_h = clm(lead_hsd, lead_couple_hsd, 0, pk)    
            self.S[pk].append_ex_mat(diag_h, upc_h, dwnc_h, ex_index)
            for s in range(self.ns):
                if tp is not None:
                    tp.log('append_lead_as_buffer(), s : {0}'.format(s))
                    tp.log('    clm()')
                diag_h, upc_h, dwnc_h = clm(lead_hsd,
                                                  lead_couple_hsd, s, pk, 'H')              
                if tp is not None:
                    tp.log('    append_ex_mat()')
                self.H[s][pk].append_ex_mat(diag_h, upc_h, dwnc_h, ex_index)                    
                if tp is not None:
                    tp.log('    clm()')
                diag_h, upc_h, dwnc_h = clm(lead_hsd,
                                                  lead_couple_hsd, s, pk, 'D')
                if tp is not None:
                    tp.log('    append_ex_mat()')
                self.D[s][pk].append_ex_mat(diag_h, upc_h, dwnc_h, ex_index, tp=tp)                 
  
    def calculate_eq_green_function(self, zp, sigma, ex=True, full=False):
        s, pk = self.s, self.pk
        #print('calculate_eq_green_function() s: {0} pk : {1}'.format(s,pk))
        self.G.reset_from_others(self.S[pk], self.H[s][pk], zp, -1, init=True)
        self.G.substract_sigma(sigma)
        if full:
            return np.linalg.inv(self.G.recover())
        else:
            #self.G.test_inv_eq()
            self.G.inv_eq()
            return self.G.recover(ex)

    def calculate_ne_green_function(self, zp, sigma, ffocc, ffvir, ex=True):
        s, pk = self.s, self.pk        
        self.G.reset_from_others(self.S[pk], self.H[s][pk], zp, -1)
        self.G.substract_sigma(sigma)
        gammaocc = []
        gammavir = []
        for ff0, ff1, tgt in zip(ffocc, ffvir, sigma):
            full_tgt = tgt.recover()
            gammaocc.append(ff0 * 1.j * (full_tgt - full_tgt.T.conj()))
            gammavir.append(ff1 * 1.j * (full_tgt - full_tgt.T.conj()))            
        glesser, ggreater = self.G.calculate_non_equilibrium_green(gammaocc,
                                                                gammavir, ex)
        return glesser, ggreater

    def abstract_sub_green_matrix(self, zp, sigma, l1, l2, inv_mat=None):
        if inv_mat is None:
            s, pk = self.s, self.pk        
            self.G.reset_from_others(self.S[pk], self.H[s][pk], zp, -1)
            self.G.substract_sigma(sigma)            
            inv_mat = self.G.inv_ne()
            gr_sub = inv_mat[l2][l1][-1]
            return gr_sub, inv_mat
        else:
            gr_sub = inv_mat[l2][l1][-1]            
            return gr_sub
       
class Tp_Sparse_Matrix:
    def __init__(self, dtype, ll_index, mat=None, band_indices=None, ex=True):
    # ll_index : lead_layer_index
    # matrix stored here will be changed to inversion

        self.lead_num = len(ll_index)
        self.ll_index = ll_index
        self.ex_ll_index = copy.deepcopy(ll_index[:])
        self.extended = ex
        self.dtype = dtype
        self.initialize()
        self.band_indices = band_indices
        if self.band_indices is None:
            self.initialize_band_indices()
        if mat is not None:
            self.reset(mat, True)
        
    def initialize_band_indices(self):
        self.band_indices = [None]
        for i in range(self.lead_num):
            self.band_indices.append([])
            for j in range(self.ex_lead_nlayer[i] - 1):
                self.band_indices[i + 1].append(None)
        
    def initialize(self):
    # diag_h : diagonal lead_hamiltonian
    # upc_h : superdiagonal lead hamiltonian
    # dwnc_h : subdiagonal lead hamiltonian 
        self.diag_h = []
        self.upc_h = []
        self.dwnc_h = []
        self.lead_nlayer = []
        self.ex_lead_nlayer = []
        self.mol_index = self.ll_index[0][0]
        self.nl = 1
        self.nb = len(self.mol_index)
        self.length = self.nb * self.nb
        self.mol_h = []

        for i in range(self.lead_num):
            self.diag_h.append([])
            self.upc_h.append([])
            self.dwnc_h.append([])
            self.lead_nlayer.append(len(self.ll_index[i]))
            if self.extended:
                self.ex_lead_nlayer.append(len(self.ll_index[i]) + 1)
            else:
                self.ex_lead_nlayer.append(len(self.ll_index[i]))
            
            assert (self.ll_index[i][0] == self.mol_index).all()
            self.nl += self.lead_nlayer[i] - 1       
            
            for j in range(self.lead_nlayer[i] - 1):
                self.diag_h[i].append([])
                self.upc_h[i].append([])
                self.dwnc_h[i].append([])
                len1 = len(self.ll_index[i][j])
                len2 = len(self.ll_index[i][j + 1])
                self.length += 2 * len1 * len2 + len2 * len2
                self.nb += len2
            
            if self.extended:                
                self.diag_h[i].append([])
                self.upc_h[i].append([])
                self.dwnc_h[i].append([])
        self.ex_nb = self.nb

    def append_ex_mat(self, diag_h, upc_h, dwnc_h, ex_index, tp=None):
        assert self.extended
        if tp is not None:
            tp.log('        append_ex_mat()')
        for i in range(self.lead_num):
            if tp is not None:
                tp.log('        append_ex_mat() i: {0}'.format(i))
            self.diag_h[i][-1] = diag_h[i]
            if tp is not None:
                tp.log('        append_ex_mat() diag_h')
            self.upc_h[i][-1] = upc_h[i]
            if tp is not None:
                tp.log('        append_ex_mat() upc_h')
            self.dwnc_h[i][-1] = dwnc_h[i]
            if tp is not None:
                tp.log('        append_ex_mat() dwnc_h')
            self.ex_ll_index[i].append(ex_index[i])
            if tp is not None:
                tp.log('        append_ex_mat() .append')
            self.ex_nb += len(ex_index[i])
            if tp is not None:
                tp.log('        append_ex_mat() += len()')
            
  
    def abstract_layer_info(self):
        self.basis_to_layer = np.empty([self.nb], int)
        self.neighbour_layers = np.zeros([self.nl, self.lead_num], int) - 1

        for i in self.mol_index:
            self.basis_to_layer[i] = 0
        nl = 1
        
        for i in range(self.lead_num):
            for j in range(self.lead_nlayer[i] - 1):
                for k in self.ll_index[i][j]:
                    self.basis_to_layer[k] = nl
                nl += 1

        nl = 1                 
        for i in range(self.lead_num):        
            self.neighbour_layers[0][i] = nl
            first = nl
            for j in range(self.lead_nlayer[i] - 1):
                if nl == first:
                    self.neighbour_layers[nl][0] = 0
                    if j != self.lead_nlayer[i] - 2:
                        self.neighbour_layers[nl][1] = nl + 1
                else:
                    self.neighbour_layers[nl][0] = nl - 1
                    if j != self.lead_nlayer[i] - 2:
                        self.neighbour_layers[nl][1] = nl + 1                    
                nl += 1
              
    def reset(self, mat, init=False):
        assert mat.dtype == self.dtype
        ind = get_matrix_index(self.mol_index)
        if init:
            self.mol_h = Banded_Sparse_Matrix(self.dtype, mat[ind.T, ind],
                                               self.band_indices[0])
            if self.band_indices[0] is None:
                self.band_indices[0] = self.mol_h.band_index            
        else:
            self.mol_h.reset(mat[ind.T, ind])

        for i in range(self.lead_num):
            for j in range(self.lead_nlayer[i] - 1):
                ind = get_matrix_index(self.ll_index[i][j])
                ind1 = get_matrix_index(self.ll_index[i][j + 1])
                indr1, indc1 = get_matrix_index(self.ll_index[i][j],
                                                      self.ll_index[i][j + 1])
                indr2, indc2 = get_matrix_index(self.ll_index[i][j + 1],
                                                  self.ll_index[i][j])
                if init:
                    self.diag_h[i][j] = Banded_Sparse_Matrix(self.dtype,
                                                             mat[ind1.T, ind1],
                                                 self.band_indices[i + 1][j])
                    if self.band_indices[i + 1][j] is None:
                        self.band_indices[i + 1][j] = \
                                                  self.diag_h[i][j].band_index
                else:
                    self.diag_h[i][j].reset(mat[ind1.T, ind1])
                self.upc_h[i][j] = mat[indr1, indc1]
                self.dwnc_h[i][j] = mat[indr2, indc2]
       
    def reset_from_others(self, tps_mm1, tps_mm2, c1, c2, init=False):
        #self.mol_h = c1 * tps_mm1.mol_h + c2 * tps_mm2.mol_h
        #print('reset_from_others {0}  {1}'.format(tps_mm1, tps_mm2))
        if init:
            self.mol_h = Banded_Sparse_Matrix(complex)
        
        self.mol_h.spar = c1 * tps_mm1.mol_h.spar + c2 * tps_mm2.mol_h.spar
        self.mol_h.band_index = tps_mm1.mol_h.band_index
        self.ex_lead_nlayer = tps_mm1.ex_lead_nlayer
        self.ex_ll_index = tps_mm1.ex_ll_index
        self.ex_nb = tps_mm1.ex_nb
        
        for i in range(self.lead_num):
            for j in range(self.ex_lead_nlayer[i] - 1):
                assert (tps_mm1.ex_ll_index[i][j] == tps_mm2.ex_ll_index[i][j]).all()
                if init:
                    self.diag_h[i][j] = Banded_Sparse_Matrix(complex)
                    self.diag_h[i][j].band_index = \
                                             tps_mm1.diag_h[i][j].band_index
                
                self.diag_h[i][j].spar = c1 * tps_mm1.diag_h[i][j].spar + \
                                      c2 * tps_mm2.diag_h[i][j].spar
                self.upc_h[i][j] = c1 * tps_mm1.upc_h[i][j] + \
                                      c2 * tps_mm2.upc_h[i][j]
                self.dwnc_h[i][j] = c1 * tps_mm1.dwnc_h[i][j] + \
                                      c2 * tps_mm2.dwnc_h[i][j]
  
    def substract_sigma(self, sigma):
        if self.extended:
            n = -2
        else:
            n = -1
        for i in range(self.lead_num):
            self.diag_h[i][n].reset_minus(sigma[i])
        
    def recover(self, ex=False):
        if ex:
            nb = self.ex_nb
            lead_nlayer = self.ex_lead_nlayer
            ll_index = self.ex_ll_index
        else:
            nb = self.nb
            lead_nlayer = self.lead_nlayer
            ll_index = self.ll_index            
        
        mat = np.zeros([nb, nb], self.dtype)
        ind = get_matrix_index(ll_index[0][0])
        
        mat[ind.T, ind] = self.mol_h.recover()
        
        gmi = get_matrix_index
        for i in range(self.lead_num):
            for j in range(lead_nlayer[i] - 1):
                ind = gmi(ll_index[i][j])
                ind1 = gmi(ll_index[i][j + 1])
                indr1, indc1 = gmi(ll_index[i][j], ll_index[i][j + 1])
                indr2, indc2 = gmi(ll_index[i][j + 1], ll_index[i][j])                
                mat[ind1.T, ind1] = self.diag_h[i][j].recover()
                mat[indr1, indc1] = self.upc_h[i][j]
                mat[indr2, indc2] = self.dwnc_h[i][j]
        return mat        

    def test_inv_eq(self, tol=1e-9):
        tp_mat = copy.deepcopy(self)
        tp_mat.inv_eq()
        mol_h = dot(tp_mat.mol_h.recover(), self.mol_h.recover())
        for i in range(self.lead_num):
            mol_h += dot(tp_mat.upc_h[i][0], self.dwnc_h[i][0])
        diff = np.max(abs(mol_h - np.eye(mol_h.shape[0])))
        if diff > tol:
            print('warning, mol_diff', diff)
        for i in range(self.lead_num):
            for j in range(self.lead_nlayer[i] - 2):
                diag_h = dot(tp_mat.diag_h[i][j].recover(),
                                                  self.diag_h[i][j].recover())
                diag_h += dot(tp_mat.dwn_h[i][j], self.upc_h[i][j])
                diag_h += dot(tp_mat.upc_h[i][j + 1], self.dwnc_h[i][j + 1])                
                diff = np.max(abs(diag_h - np.eye(diag_h.shape[0])))
                if diff > tol:
                    print('warning, diag_diff', i, j, diff)
            j = self.lead_nlayer[i] - 2
            diag_h = dot(tp_mat.diag_h[i][j].recover(),
                                                  self.diag_h[i][j].recover())
            diag_h += dot(tp_mat.dwnc_h[i][j], self.upc_h[i][j])
            diff = np.max(abs(diag_h - np.eye(diag_h.shape[0])))
            if diff > tol:
                print('warning, diag_diff', i, j, diff)            
                                                
    def inv_eq(self):
        q_mat = []
        for i in range(self.lead_num):
            q_mat.append([])
            nll = self.lead_nlayer[i]
            for j in range(nll - 1):
                q_mat[i].append([])
            end = nll - 2
            q_mat[i][end] =  self.diag_h[i][end].inv()
          
            for j in range(end - 1, -1, -1):
                self.diag_h[i][j].reset_minus(self.dotdot(
                                                    self.upc_h[i][j + 1],
                                                         q_mat[i][j + 1],
                                            self.dwnc_h[i][j + 1]), full=True)
                q_mat[i][j] = self.diag_h[i][j].inv()
        h_mm = self.mol_h

        for i in range(self.lead_num):
            h_mm.reset_minus(self.dotdot(self.upc_h[i][0], q_mat[i][0],
                                                self.dwnc_h[i][0]), full=True)
        inv_h_mm = h_mm.inv()
        h_mm.reset(inv_h_mm)
        
        for i in range(self.lead_num):
            tmp_dc = self.dwnc_h[i][0].copy()
            self.dwnc_h[i][0] = -self.dotdot(q_mat[i][0], tmp_dc, inv_h_mm)
            self.upc_h[i][0] = -self.dotdot(inv_h_mm, self.upc_h[i][0],
                                                                    q_mat[i][0])
            dim = len(self.ll_index[i][1])
            self.diag_h[i][0].reset(dot(q_mat[i][0], np.eye(dim) -
                                                  dot(tmp_dc, self.upc_h[i][0])))
            for j in range(1, self.lead_nlayer[i] - 1):
                tmp_dc = self.dwnc_h[i][j].copy()
                self.dwnc_h[i][j] = -self.dotdot(q_mat[i][j], tmp_dc,
                                                self.diag_h[i][j - 1].recover())
                self.upc_h[i][j] = -self.dotdot(self.diag_h[i][j - 1].recover(),
                                                    self.upc_h[i][j],
                                                     q_mat[i][j])
                dim = len(self.ll_index[i][j + 1])
                self.diag_h[i][j].reset(dot(q_mat[i][j], np.eye(dim) -
                                           dot(tmp_dc, self.upc_h[i][j])))

    def inv_ne(self):
        q_mat = []
        qi_mat = []
        inv_mat = []
        #structure of inv_mat inv_cols_1, inv_cols_2, ..., inv_cols_n (n:lead_num)
        #structure of inv_cols_i   inv_cols_l1, inv_cols_l2,..., inv_cols_ln, inv_cols_mm(matrix)
        #structure of inv_cols_li  inv_cols_ll1, inv_cols_ll2,...,inv_cols_ll3
        for i in range(self.lead_num):
            q_mat.append([])
            qi_mat.append([])
            inv_mat.append([])
            
            nll = self.lead_nlayer[i]
            for j in range(nll - 1):
                q_mat[i].append([])
                qi_mat[i].append([])
                
            for j in range(self.lead_num):
                inv_mat[i].append([])
                nll_j = self.lead_nlayer[j]
                for k in range(nll_j - 1):
                    inv_mat[i][j].append([])
            inv_mat[i].append([])                
            
            end = nll - 2
            q_mat[i][end] =  self.diag_h[i][end].inv()
            for j in range(end - 1, -1, -1):
                tmp_diag_h = copy.deepcopy(self.diag_h[i][j])
                tmp_diag_h.reset_minus(self.dotdot(self.upc_h[i][j + 1],
                                                     q_mat[i][j + 1],
                                                  self.dwnc_h[i][j + 1]),
                                        full=True)
                q_mat[i][j] = tmp_diag_h.inv()
        # above get all the q matrix, then if want to solve the cols
        # cooresponding to the lead i, the q_mat[i] will not be used

        #q_mm = self.mol_h.recover()
        q_mm = copy.deepcopy(self.mol_h)
        for i in range(self.lead_num):
            #q_mm -= dot(dot(self.upc_h[i][0], q_mat[i][0]),
            #                                 self.dwnc_h[i][0])
            q_mm.reset_minus(self.dotdot(self.upc_h[i][0],
                                  q_mat[i][0], self.dwnc_h[i][0]), full=True)
        
        for i in range(self.lead_num):
        # solve the corresponding cols to the lead i
            nll = self.lead_nlayer[i]
    
            #qi_mat[i][0] = q_mm + self.dotdot(self.upc_h[i][0],q_mat[i][0],
            #                                                self.dwnc_h[i][0])
            q_mm_tmp = copy.deepcopy(q_mm)
            q_mm_tmp.reset_plus(self.dotdot(self.upc_h[i][0],q_mat[i][0],
                                                self.dwnc_h[i][0]), full=True)
            
            #inv(qi_mat[i][0])
            qi_mat[i][0] = q_mm_tmp.inv()
            for j in range(1, nll - 1):
                tmp_diag_h = copy.deepcopy(self.diag_h[i][j - 1])
                tmp_diag_h.reset_minus(self.dotdot(self.dwnc_h[i][j -1],
                                                        qi_mat[i][j - 1],
                                                        self.upc_h[i][j - 1]),
                                                       full=True)
                qi_mat[i][j] = tmp_diag_h.inv()

            tmp_diag_h = copy.deepcopy(self.diag_h[i][nll - 2])
            tmp_diag_h.reset_minus(self.dotdot(self.dwnc_h[i][nll - 2],
                                                qi_mat[i][nll -2],
                                               self.upc_h[i][nll -2]),
                                                  full=True)
            inv_mat[i][i][nll - 2] = tmp_diag_h.inv()
            
            for j in range(nll - 3, -1, -1):
                inv_mat[i][i][j] = -self.dotdot(qi_mat[i][j + 1],
                                                  self.upc_h[i][j + 1],
                                                 inv_mat[i][i][j + 1])
            inv_mat[i][self.lead_num] = -self.dotdot(qi_mat[i][0],
                                                  self.upc_h[i][0],
                                                  inv_mat[i][i][0]) 
            
            for j in range(self.lead_num):
                if j != i:
                    nlj = self.lead_nlayer[j]
                    inv_mat[i][j][0] = -self.dotdot(q_mat[j][0], self.dwnc_h[j][0],
                                                inv_mat[i][self.lead_num])
                    for k in range(1, nlj - 1):
                        inv_mat[i][j][k] = -self.dotdot(q_mat[j][k], self.dwnc_h[j][k],
                                                inv_mat[i][j][k - 1])                         
        return inv_mat 
  
  
    def combine_inv_mat(self, inv_mat):
        nb = self.nb
        mat = np.zeros([nb, nb], complex)
        for i in range(self.lead_num):
            indr, indc = get_matrix_index(self.ll_index[i][0],
                                          self.ll_index[i][-1])
            mat[indr, indc] = inv_mat[i][self.lead_num]
            for j in range(self.lead_num):
                for k in range(1, self.lead_nlayer[j]):
                    indr, indc = get_matrix_index(self.ll_index[j][k],
                                                  self.ll_index[i][-1])
                    mat[indr, indc] = inv_mat[i][j][k - 1]
        return mat
  
    def dotdot(self, mat1, mat2, mat3):
        return dot(mat1, dot(mat2, mat3))
    
    def calculate_non_equilibrium_green(self, se_less, se_great, ex=True):
        inv_mat = self.inv_ne()
        glesser = self.calculate_keldysh_green(inv_mat, se_less, ex)
        ggreater = self.calculate_keldysh_green(inv_mat, se_great, ex)
        return glesser, ggreater  
    
    def calculate_keldysh_green(self, inv_mat, keldysh_se, ex=True):
        #se_less less selfenergy, structure  se_1, se_2, se_3,..., se_n
        #the lead sequence of se_less should be the same to self.ll_index
        self.mol_h.spar.fill(0.0)
        for i in range(self.lead_num):
            nll = self.lead_nlayer[i]
            for j in range(nll - 1):
                self.diag_h[i][j].spar.fill(0.0)
                self.upc_h[i][j].fill(0.0)
                self.dwnc_h[i][j].fill(0.0)                    
        
        for i in range(self.lead_num):
            # less selfenergy loop
            self.mol_h.reset_plus(self.dotdot(inv_mat[i][self.lead_num],
                                              keldysh_se[i],
                                          inv_mat[i][self.lead_num].T.conj()),
                                              full=True)            
            for j in range(self.lead_num):
               # matrix operation loop    
                nlj = self.lead_nlayer[j]
                self.diag_h[j][0].reset_plus(self.dotdot(inv_mat[i][j][0],
                                                    keldysh_se[i],
                                                 inv_mat[i][j][0].T.conj()),
                                                      full=True)            
            
                self.dwnc_h[j][0] += self.dotdot(inv_mat[i][j][0], keldysh_se[i],
                                           inv_mat[i][self.lead_num].T.conj())
            
                self.upc_h[j][0] += self.dotdot(inv_mat[i][self.lead_num],
                                                keldysh_se[i],
                                            inv_mat[i][j][0].T.conj())
            
                for k in range(1, nlj -1):
                    self.diag_h[j][k].reset_plus(self.dotdot(inv_mat[i][j][k],
                                                             keldysh_se[i],
                                                 inv_mat[i][j][k].T.conj()),
                                                    full=True)
                    
                    self.dwnc_h[j][k] += self.dotdot(inv_mat[i][j][k],
                                                     keldysh_se[i],
                                                inv_mat[i][j][k - 1].T.conj())
                        
                    self.upc_h[j][k] +=  self.dotdot(inv_mat[i][j][k - 1],
                                                     keldysh_se[i],
                                                    inv_mat[i][j][k].T.conj())
        return self.recover(ex)            

    def test_inv_speed(self):
        full_mat = self.recover()
        timer = Timer()
        timer.start('full_numpy')
        tmp0 = np.linalg.inv(full_mat)
        timer.stop('full_numpy')
        
        timer.start('full_lapack')
        inverse_general(full_mat)
        timer.stop('full_lapack')
        
        timer.start('sparse_lapack')
        self.inv_eq()
        timer.stop('sparse_lapack')
        
        timer.start('sparse_lapack_ne')
        self.inv_ne()
        timer.stop('sparse_lapack_ne')
        
        times = []
        methods = ['full_numpy', 'full_lapack', 'sparse_lapack']
        for name in methods:
            time = timer.timers[name,]
            print(name, time)
            times.append(time)
        
        mintime = np.min(times)
        self.inv_method = methods[np.argmin(times)]
        print('mintime', mintime)
        
        print('sparse_lapack_ne', timer.timers['sparse_lapack_ne',])

class CP_Sparse_HSD:
    def __init__(self, dtype, ns, npk, index=None):
        self.index = index
        self.dtype = dtype
        self.H = []
        self.S = []
        self.D = []
        self.ns = ns
        self.npk = npk
        self.s = 0
        self.pk = 0
        for s in range(ns):
            self.H.append([])
            self.D.append([])
            for k in range(npk):
                self.H[s].append([])
                self.D[s].append([])
        for k in range(npk):
            self.S.append([])        

    def reset(self, s, pk, mat, flag='S', init=False):
        assert mat.dtype == self.dtype
        if flag == 'S':
            spar = self.S
        elif flag == 'H':
            spar = self.H[s]
        elif flag == 'D':
            spar = self.D[s]
        if not init:
            spar[pk].reset(mat)
        elif self.index is not None:
            spar[pk] = CP_Sparse_Matrix(self.dtype, mat, self.index)
        else:
            spar[pk] = CP_Sparse_Matrix(self.dtype, mat)
            self.index = spar[pk].index
     
class CP_Sparse_Matrix:
    def __init__(self, dtype, mat=None, index=None, flag=None, tol=1e-9):
        self.tol = tol
        self.index = index
        self.dtype = dtype
        self.flag = flag
        if mat is not None:
            if self.index is None:
                self.initialize(mat)
            else:
                self.reset(mat)
        
    def initialize(self, mat):
        assert self.dtype == mat.dtype
        dim = mat.shape[-1]
        ud_array = np.empty([dim])
        dd_array = np.empty([dim])
        for i in range(dim):
            ud_array[i] = np.sum(abs(np.diag(mat, i)))
            dd_array[i] = np.sum(abs(np.diag(mat, -i)))
        spar_sum = 0
        mat_sum = np.sum(abs(mat))
        if np.sum(abs(ud_array)) >  np.sum(abs(dd_array)):
            self.flag = 'U'
            i = -1
            while abs(mat_sum - spar_sum) > self.tol * 10:
                i += 1
                while ud_array[i] < self.tol and  i < dim - 1:
                    i += 1
                self.index = (i, dim)
                ldab = dim - i
                self.spar = mat[:ldab, i:].copy()
                spar_sum = np.sum(abs(self.spar))
        else:
            self.flag = 'L'
            i = -1
            while abs(mat_sum - spar_sum) > self.tol * 10:
                i += 1
                while dd_array[i] < self.tol and  i < dim - 1:
                    i += 1
                self.index = (-i, dim)
                ldab = dim - i
                self.spar = mat[i:, :ldab].copy()
                spar_sum = np.sum(abs(self.spar))
        
    def reset(self, mat):
        assert mat.dtype == self.dtype and mat.shape[-1] == self.index[1]
        dim = mat.shape[-1]
        if self.index[0] > 0:
            ldab = dim - self.index[0]
            self.spar = mat[:ldab, self.index[0]:].copy()            
        else:
            ldab = dim + self.index[0]
            self.spar = mat[-self.index[0]:, :ldab].copy()               

    def recover(self, trans='n'):
        nb = self.index[1]
        mat = np.zeros([nb, nb], self.dtype)
        if self.index[0] > 0:
            ldab = nb - self.index[0]
            mat[:ldab, self.index[0]:] = self.spar
        else:
            ldab = nb + self.index[0]
            mat[-self.index[0]:, :ldab] = self.spar
        if trans == 'c':
            if self.dtype == float:
                mat = mat.T.copy()
            else:
                mat = mat.T.conj()
        return mat

class Se_Sparse_Matrix:
    def __init__(self, mat, tri_type, nn=None, tol=1e-9):
        # coupling sparse matrix A_ij!=0 if i>dim-nn and j>dim-nn (for right selfenergy)
        # or A_ij!=0 if i<nn and j<nn (for left selfenergy, dim is the shape of A)
        self.tri_type = tri_type
        self.tol = tol
        self.nb = mat.shape[-1]
        self.spar = []
        if nn is None:
            self.initialize(mat)
        else:
            self.reset(mat, nn)

    def initialize(self, mat):
        self.nn = 0
        nb = self.nb
        tol = self.tol
        if self.tri_type == 'L':
            while self.nn < nb and np.sum(abs(mat[self.nn])) > tol:
                self.nn += 1
            self.spar = mat[:self.nn, :self.nn].copy()
        else:
            while self.nn < nb and np.sum(abs(mat[nb - self.nn - 1])) > tol:
                self.nn += 1
            self.spar = mat[-self.nn:, -self.nn:].copy()                 
        
        diff = abs(np.sum(abs(mat)) - np.sum(abs(self.spar)))
        if diff > tol * 10:
            print('Warning! Sparse Matrix Diff', diff)
        
    def reset(self, mat, nn=None):
        if nn is not None:
            self.nn = nn
        if self.tri_type == 'L':
            self.spar = mat[:self.nn, :self.nn].copy()
        else:
            self.spar = mat[-self.nn:, -self.nn:].copy()    
  
    def restore(self):
        mat = np.zeros([self.nb, self.nb], complex)
        if self.tri_type == 'L':
            mat[:self.nn, :self.nn] = self.spar
        else:
            mat[-self.nn:, -self.nn:] = self.spar
        return mat   

