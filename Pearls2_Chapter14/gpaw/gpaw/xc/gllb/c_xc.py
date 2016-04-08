from gpaw.xc.gllb.contribution import Contribution
from gpaw.xc import XC
from gpaw.xc.pawcorrection import rnablaY_nLv
from gpaw.sphere.lebedev import weight_n
import numpy as np
from numpy import dot as dot3  # Avoid dotblas bug!
from math import pi, sqrt

class C_XC(Contribution):
    def __init__(self, nlfunc, weight, functional = 'LDA'):
        Contribution.__init__(self, nlfunc, weight)
        self.functional = functional

    def get_name(self):
        return 'XC'

    def get_desc(self):
        return "("+self.functional+")"
        
    def initialize(self):
        self.xc = XC(self.functional)
        self.vt_sg = self.nlfunc.finegd.empty(self.nlfunc.nspins)
        self.e_g = self.nlfunc.finegd.empty()

    def initialize_1d(self):
        self.ae = self.nlfunc.ae
        self.xc = XC(self.functional) 
        self.v_g = np.zeros(self.ae.N)

    def calculate_spinpaired(self, e_g, n_g, v_g):
        self.e_g[:] = 0.0
        self.vt_sg[:] = 0.0
        self.xc.calculate(self.nlfunc.finegd, n_g[None, ...], self.vt_sg,
                          self.e_g)
        v_g += self.weight * self.vt_sg[0]
        e_g += self.weight * self.e_g

    def calculate_spinpolarized(self, e_g, n_sg, v_sg):
        self.e_g[:] = 0.0
        self.vt_sg[:] = 0.0
        self.xc.calculate(self.nlfunc.finegd, n_sg, self.vt_sg, self.e_g)
        #self.xc.get_energy_and_potential(na_g, self.vt_sg[0], nb_g, self.vt_sg[1], e_g=self.e_g)
        v_sg[0] += self.weight * self.vt_sg[0]
        v_sg[1] += self.weight * self.vt_sg[1]
        e_g += self.weight * self.e_g

    def calculate_energy_and_derivatives(self, setup, D_sp, H_sp, a, addcoredensity=True):
        E = self.xc.calculate_paw_correction(setup, D_sp, H_sp, True, a)
        E += setup.xc_correction.Exc0
        print("E", E)
        return E

    def add_xc_potential_and_energy_1d(self, v_g):
        self.v_g[:] = 0.0
        Exc = self.xc.calculate_spherical(self.ae.rgd,
                                          self.ae.n.reshape((1, -1)),
                                          self.v_g.reshape((1, -1)))
        v_g += self.weight * self.v_g
        return self.weight * Exc

    def add_smooth_xc_potential_and_energy_1d(self, vt_g):
        self.v_g[:] = 0.0
        Exc = self.xc.calculate_spherical(self.ae.rgd,
                                          self.ae.nt.reshape((1, -1)),
                                          self.v_g.reshape((1, -1)))
        vt_g += self.weight * self.v_g
        return self.weight * Exc

    def initialize_from_atomic_orbitals(self, basis_functions):
        # LDA needs only density, which is already initialized
        pass

    def add_extra_setup_data(self, dict):
        # LDA has not any special data
        pass

    def write(self, writer, natoms):
        # LDA has not any special data to be written
        pass

    def read(self, reader):
        # LDA has not any special data to be read
        pass
        
