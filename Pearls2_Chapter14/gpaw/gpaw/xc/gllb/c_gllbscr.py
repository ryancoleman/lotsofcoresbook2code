from gpaw.xc.gllb.contribution import Contribution
from gpaw.xc import XC
from gpaw.xc.pawcorrection import rnablaY_nLv
from gpaw.xc.gllb import safe_sqr
from math import sqrt, pi, exp
from gpaw.utilities import erf
from gpaw.io.tar import TarFileReference
from gpaw.sphere.lebedev import weight_n
import numpy as np

K_G = 0.382106112167171

class C_GLLBScr(Contribution):
    def __init__(self, nlfunc, weight, functional='GGA_X_B88', width=None, eps=0.05, damp=1e-10):
        Contribution.__init__(self, nlfunc, weight)
        self.functional = functional
        self.old_coeffs = None
        self.iter = 0
        self.damp = damp
        if width is not None:
            width = width / 27.21
        self.eps = eps / 27.21
        self.width = width
 
    def get_name(self):
        return 'SCREENING'

    def set_damp(self, damp):
        self.damp = damp

    def get_desc(self):
        return '(' + self.functional + ')'
        
    # Initialize GLLBScr functional
    def initialize_1d(self):
        self.ae = self.nlfunc.ae
        self.xc = XC(self.functional)
        self.v_g = np.zeros(self.ae.N)
        self.e_g = np.zeros(self.ae.N)

    # Calcualte the GLLB potential and energy 1d
    def add_xc_potential_and_energy_1d(self, v_g):
        self.v_g[:] = 0.0
        self.e_g[:] = 0.0
        self.xc.calculate_spherical(self.ae.rgd, self.ae.n.reshape((1, -1)),
                                    self.v_g.reshape((1, -1)), self.e_g)
        v_g += 2 * self.weight * self.e_g / (self.ae.n + self.damp)
        Exc = self.weight * np.sum(self.e_g * self.ae.rgd.dv_g)
        return Exc

    def initialize(self):
        self.occupations = self.nlfunc.occupations
        self.xc = XC(self.functional)

        # Always 1 spin, no matter what calculation nspins is
        self.vt_sg = self.nlfunc.finegd.empty(1) 
        self.e_g = self.nlfunc.finegd.empty()#.ravel()

    def get_coefficient_calculator(self):
        return self

    def f(self, f):
        if self.width is None:
            if f > self.eps:
                return sqrt(f)
            else:
                return 0.0
        else:
            dEH = -f
            w = self.width
            if dEH / w < -100:
                return sqrt(f)
            Knew = -0.5 * erf(sqrt((max(0.0,dEH)-dEH)/w)) * \
                    sqrt(w*pi) * exp(-dEH/w)
            Knew += 0.5 * sqrt(w*pi)*exp(-dEH/w)
            Knew += sqrt(max(0.0,dEH)-dEH)*exp(max(0.0,dEH)/w)
            #print dEH, w, dEH/w, Knew, f**0.5
            return Knew
    
    def get_coefficients(self, e_j, f_j):
        homo_e = max( [ np.where(f>1e-3, e, -1000) for f,e in zip(f_j, e_j)] ) 
        return [ f * K_G * self.f(homo_e - e) for e,f in zip(e_j, f_j) ]

    def get_coefficients_1d(self, smooth=False, lumo_perturbation = False):
        homo_e = max( [ np.where(f>1e-3, e, -1000) for f,e in zip(self.ae.f_j, self.ae.e_j)]) 
        if not smooth:
            if lumo_perturbation:
                lumo_e = min( [ np.where(f<1e-3, e, 1000) for f,e in zip(self.ae.f_j, self.ae.e_j)])
                return np.array([ f * K_G * (self.f(lumo_e - e) - self.f(homo_e -e))
                                        for e,f in zip(self.ae.e_j, self.ae.f_j) ])
            else:
                return np.array([ f * K_G * (self.f(homo_e - e))
                                   for e,f in zip(self.ae.e_j, self.ae.f_j) ])
        else:
            return [ [ f * K_G * self.f(homo_e - e)
                    for e,f in zip(e_n, f_n) ]
                     for e_n, f_n in zip(self.ae.e_ln, self.ae.f_ln) ]
        

    def get_coefficients_by_kpt(self, kpt_u, lumo_perturbation=False, homolumo=None, nspins=1):
        #if not hasattr(kpt_u[0],'orbitals_ready'):
        #    kpt_u[0].orbitals_ready = True
        #    return None
        #if not hasattr(self.occupations, 'nvalence'):
        #    print "occupations not ready"
        #    return None
        #if self.occupations.nvalence is None:
        #    return None
        #if kpt_u[0].psit_nG is None or isinstance(kpt_u[0].psit_nG,
        #                                          TarFileReference): 
        #    if kpt_u[0].C_nM is None:
        #        return None

        if homolumo is None:
            # Find homo and lumo levels for each spin
            eref_s = []
            eref_lumo_s = []  
            for s in range(nspins):
                homo, lumo = self.occupations.get_homo_lumo_by_spin(self.nlfunc.wfs, s)
                eref_s.append(homo)
                eref_lumo_s.append(lumo)
        else:
            eref_s, eref_lumo_s = homolumo
            if not isinstance(eref_s, (list, tuple)):
                eref_s = [ eref_s ]
                eref_lumo_s = [ eref_lumo_s ]

        # The parameter ee might sometimes be set to small thereshold value to
        # achieve convergence on small systems with degenerate HOMO.
        if len(kpt_u) > nspins:
            ee = 0.0
        else:
            ee = 0.05 / 27.21

        if lumo_perturbation:
            return [np.array([
                f * K_G * (self.f(eref_lumo_s[kpt.s]-e)
                          -self.f(eref_s[kpt.s]-e))
                     for e, f in zip(kpt.eps_n, kpt.f_n) ])
                     for kpt in kpt_u ]
            
            
        else:
            coeff = [ np.array([ f * K_G * self.f(eref_s[kpt.s] - e) 
                     for e, f in zip(kpt.eps_n, kpt.f_n) ])
                     for kpt in kpt_u ]
            #print coeff
            return coeff
        

    def calculate_spinpaired(self, e_g, n_g, v_g):
        self.e_g[:] = 0.0
        self.vt_sg[:] = 0.0
        self.xc.calculate(self.nlfunc.finegd, n_g[None, ...], self.vt_sg,
                          self.e_g)
        self.e_g[:] = np.where(n_g<self.damp, 0, self.e_g)
        v_g += self.weight * 2 * self.e_g / (n_g + self.damp)
        e_g += self.weight * self.e_g

    def calculate_spinpolarized(self, e_g, n_sg, v_sg):
        # Calculate spinpolarized exchange screening as two spin-paired calculations n=2*n_s
        for n, v in [ (n_sg[0], v_sg[0]), (n_sg[1], v_sg[1]) ]:
            self.e_g[:] = 0.0
            self.vt_sg[:] = 0.0
            self.xc.calculate(self.nlfunc.finegd, 2*n[None, ...], self.vt_sg, self.e_g)
            self.e_g[:] = np.where(n<self.damp, 0, self.e_g)
            v += self.weight * 2 * self.e_g / (2 * n + self.damp)
            e_g += self.weight * self.e_g / 2

    def calculate_energy_and_derivatives(self, setup, D_sp, H_sp, a, addcoredensity=True):
        # Get the XC-correction instance
        c = setup.xc_correction
        nspins = self.nlfunc.nspins

        E = 0
        for D_p, dEdD_p in zip(D_sp, H_sp):
                D_Lq = np.dot(c.B_pqL.T, nspins*D_p)
                n_Lg = np.dot(D_Lq, c.n_qg)
                if addcoredensity:
                     n_Lg[0] += c.nc_g * sqrt(4 * pi)
                nt_Lg = np.dot(D_Lq, c.nt_qg)
                if addcoredensity:
                     nt_Lg[0] += c.nct_g * sqrt(4 * pi)
                dndr_Lg = np.zeros((c.Lmax, c.ng))
                dntdr_Lg = np.zeros((c.Lmax, c.ng))
                for L in range(c.Lmax):
                    c.rgd.derivative(n_Lg[L], dndr_Lg[L])
                    c.rgd.derivative(nt_Lg[L], dntdr_Lg[L])
                vt_g = np.zeros(c.ng)
                v_g = np.zeros(c.ng)
                e_g = np.zeros(c.ng)
                deda2_g = np.zeros(c.ng)
                for y, (w, Y_L) in enumerate(zip(weight_n, c.Y_nL)):
                    # Cut gradient releated coefficient to match the setup's Lmax
                    A_Li = rnablaY_nLv[y, :c.Lmax]
                    
                    # Expand pseudo density
                    nt_g = np.dot(Y_L, nt_Lg)
                    
                    # Expand pseudo density gradient
                    a1x_g = np.dot(A_Li[:, 0], nt_Lg)
                    a1y_g = np.dot(A_Li[:, 1], nt_Lg)
                    a1z_g = np.dot(A_Li[:, 2], nt_Lg)
                    a2_g = a1x_g**2 + a1y_g**2 + a1z_g**2
                    a2_g[1:] /= c.rgd.r_g[1:]**2
                    a2_g[0] = a2_g[1]
                    a1_g = np.dot(Y_L, dntdr_Lg)
                    a2_g += a1_g**2
                   
                    vt_g[:] = 0.0
                    e_g[:] = 0.0
                    # Calculate pseudo GGA energy density (potential is discarded)
                    self.xc.kernel.calculate(e_g, nt_g.reshape((1, -1)),
                                             vt_g.reshape((1, -1)),
                                         a2_g.reshape((1, -1)),
                                         deda2_g.reshape((1, -1)))

                    # Calculate pseudo GLLB-potential from GGA-energy density
                    vt_g[:] = 2 * e_g / (nt_g + self.damp)

                    dEdD_p -= self.weight * w * np.dot(np.dot(c.B_pqL, Y_L),
                                          np.dot(c.nt_qg, vt_g * c.rgd.dv_g))

                    E -= w * np.dot(e_g, c.rgd.dv_g) / nspins

                    # Expand density
                    n_g = np.dot(Y_L, n_Lg)

                    # Expand density gradient
                    a1x_g = np.dot(A_Li[:, 0], n_Lg)
                    a1y_g = np.dot(A_Li[:, 1], n_Lg)
                    a1z_g = np.dot(A_Li[:, 2], n_Lg)
                    a2_g = a1x_g**2 + a1y_g**2 + a1z_g**2
                    a2_g[1:] /= c.rgd.r_g[1:]**2
                    a2_g[0] = a2_g[1]
                    a1_g = np.dot(Y_L, dndr_Lg)
                    a2_g += a1_g**2

                    v_g[:] = 0.0
                    e_g[:] = 0.0
                    # Calculate GGA energy density (potential is discarded)
                    self.xc.kernel.calculate(e_g, n_g.reshape((1, -1)),
                                             v_g.reshape((1, -1)),
                                             a2_g.reshape((1, -1)),
                                             deda2_g.reshape((1, -1)))

                    # Calculate GLLB-potential from GGA-energy density
                    v_g[:] = 2 * e_g / (n_g + self.damp)

                    dEdD_p += self.weight * w * np.dot(np.dot(c.B_pqL, Y_L),
                                          np.dot(c.n_qg, v_g * c.rgd.dv_g))
                    E += w * np.dot(e_g, c.rgd.dv_g) / nspins

        return E * self.weight

    def add_smooth_xc_potential_and_energy_1d(self, vt_g):
        self.v_g[:] = 0.0
        self.e_g[:] = 0.0
        self.xc.calculate_spherical(self.ae.rgd, self.ae.nt.reshape((1, -1)),
                                    self.v_g.reshape((1, -1)), self.e_g)
        vt_g += 2 * self.weight * self.e_g / (self.ae.nt + self.damp)
        return self.weight * np.sum(self.e_g * self.ae.rgd.dv_g)

    def initialize_from_atomic_orbitals(self, basis_functions):
        # GLLBScr needs only density which is already initialized
        pass
        
    def add_extra_setup_data(self, dict):
        # GLLBScr has not any special data
        pass

    def read(self, reader):
        # GLLBScr has no special data to be read
        pass

    def write(self, writer, natoms):
        # GLLBScr has no special data to be written
        pass
        


