from gpaw.transport.tools import dagger, dot
from gpaw.transport.sparse_matrix import Banded_Sparse_Matrix, Se_Sparse_Matrix
import copy
import numpy as np
import cPickle

class LeadSelfEnergy:
    #This object use the sparse_matrix object Banded_Sparse_HSD
    conv = 1e-8 # Convergence criteria for surface Green function
    def __init__(self, hsd_ii, hsd_ij, data_path=None, direction='left',
                 rotation_mat=None):
        self.hsd_ii = hsd_ii
        self.hsd_ij = hsd_ij
        self.data_path = data_path
        self.direction = direction
        self.energy = None
        self.rotation_mat = rotation_mat
        self.bias = 0
        self.s = 0
        self.pk = 0
        self.nid_plus = 30
        self.sigma = None
        self.tol = 1e-6
        
    def __call__(self, energy, flag=None):
        if self.data_path is not None and flag is not None:
            nid = int(flag[2:])
            nid += self.nid_plus
            flag = str(flag[:2]) + str(nid)
            fd = file(self.data_path[self.direction] + '/' +
                      self.direction + '/' + flag, 'r')
            data = cPickle.load(fd)
            return data
        else:
            if self.energy is None or abs(energy - self.energy) > self.tol:
                self.energy = energy
                z = energy - self.bias         
                tau_ij = z * self.hsd_ij.S[self.pk].recover() - \
                                     self.hsd_ij.H[self.s][self.pk].recover()
                tau_ji = z * dagger(self.hsd_ij.S[self.pk].recover()) - \
                             dagger(self.hsd_ij.H[self.s][self.pk].recover())
                ginv = self.get_sgfinv(energy)
                a_ij = dot(ginv, tau_ij)
                mat = dot(tau_ji, a_ij)
                #if self.rotation_mat is not None:
                #    mat = np.dot(self.rotation_mat, mat)
                #    mat = np.dot(mat, self.rotation_mat.T)
                self.sigma = Banded_Sparse_Matrix(complex, mat,
                                            self.hsd_ii.S[self.pk].band_index)
            return self.sigma
       
    def set_bias(self, bias):
        self.bias = bias
        self.nid_plus = -int(bias // 0.05) + 30
        
    def get_lambda(self, energy):
        sigma = self(energy)
        sigma_mm = sigma.recover()
        return 1.j * (sigma_mm - dagger(sigma_mm))
    
    def get_sgfinv(self, energy):
        """The inverse of the retarded surface Green function"""
        z = energy - self.bias 
        v_00 = Banded_Sparse_Matrix(complex,
                                     None,
                                     self.hsd_ii.S[self.pk].band_index)
        
        v_00.reset_from_others(self.hsd_ii.S[self.pk],
                                self.hsd_ii.H[self.s][self.pk],
                                z, -1.0)
        
        v_11 = copy.deepcopy(v_00)
        
        v_10 = z * self.hsd_ij.S[self.pk].recover()- \
                                     self.hsd_ij.H[self.s][self.pk].recover()
        v_01 = z * dagger(self.hsd_ij.S[self.pk].recover()) - \
                              dagger(self.hsd_ij.H[self.s][self.pk].recover())
        delta = self.conv + 1
        while delta > self.conv:
            inv_v_11 = v_11.inv()
            a = dot(inv_v_11, v_01)
            b = dot(inv_v_11, v_10)
            v_01_dot_b = dot(v_01, b)
            v_00.reset_minus(v_01_dot_b, full=True)
            v_11.reset_minus(dot(v_10, a), full=True)
            v_11.reset_minus(v_01_dot_b, full=True)
            v_01 = -dot(v_01, a)
            v_10 = -dot(v_10, b)
            delta = np.abs(v_01).max()
        return v_00.inv()


