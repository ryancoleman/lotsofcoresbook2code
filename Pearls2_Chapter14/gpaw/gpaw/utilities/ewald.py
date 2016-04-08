import numpy as np

from gpaw.utilities import erf
from math import pi


def erfc(x):
    """The complimentary error function."""
    return 1. - erf(x)


class Ewald:
    """Class for calculating Ewald summations.
    
    'cell' is the unit cell in the cartesian basis.
    
    'G' is a normalized Ewald parameter. Large G results in fast
    real-space convergence but slow convergence in reciprocal space.
    This describes the width (in reciprocal space, and inverse width
    in real space) of the Gaussian shaped probe charge.
    
    Ng and Nl are lists specifying the number or nearest neighbors in
    sums over the reciprocal lattice and the real space lattice
    respectively.
    """
    
    def __init__(self, cell, G=5, Ng=[9, 9, 9], Nl=[3, 3, 3]):
        self.cell = cell
        self.Vcell = abs(np.linalg.det(cell))
        self.recip_cell = np.linalg.inv(self.cell).T
        self.Ng = Ng
        self.Nl = Nl
        self.G = G

        # Renormalize G to the longest lattice vector
        self.G /= np.sqrt((cell**2).sum(1).max())

    def get_wigner_seitz_radius(self, N):
        """Wigner-Seitz radius for N electrons in the unit cell."""
        return ( 3. * self.Vcell / (4. * pi * N))**(1./3.)

    def get_sum_recip_ij(self, r_v, eps=1e-10):
        """The reciprocal space sum.

        ::
        
              -----  -x  i g.r
           pi \     e   e
          ---2 |   -------   , with x = g^2/(4 G^2)
          V G /       x
              -----
             g not 0
        """
        N = self.Ng
        b1, b2, b3 = self.recip_cell
        E_g = 0.
        for i in range(-N[0], N[0] + 1):
            for j in range(-N[1], N[1] + 1):
                for k in range(-N[2], N[2] + 1):
                    g_v = 2 * pi * (i * b1 + j * b2 + k * b3)
                    g2 = np.dot(g_v, g_v)
                    if g2 > eps: # exclude g=0
                        x = .25 * g2 / self.G**2
                        E_g += np.exp(1.j * np.dot(r_v, g_v) - x) / x
        return pi / (self.Vcell * self.G**2) * E_g.real   

    def get_sum_real_ij(self, r_v, eps=1e-5):
        """The real space sum.

        ::
        
            -----   
            \     erfc( G [l-r| )
             |   -------------- 
            /         |l-r|
            -----
           l not 0
           
        Note: Add the l=0 term with erfc(r).
        """
        N = self.Nl
        a1, a2, a3 = self.cell
        E_r = 0.
        for i in range(-N[0], N[0] + 1):
            for j in range(-N[1], N[1] + 1):
                for k in range(-N[2], N[2] + 1):       
                    l_v = i * a1 + j * a2 + k * a3
                    if np.linalg.norm(l_v) > eps: # exclude l=0
                        lr = np.linalg.norm(l_v - r_v)
                        E_r += erfc(self.G * lr) / lr
        return E_r

    def get_electrostatic_potential(self, r, r_B, q_B, excludefroml0=None):
        """...
        
        Calculates the electrostatic potential at point r_i from point
        charges at {r_B} in a lattice using the Ewald summation.
        
        Charge neutrality is obtained by adding the homogenous density
        q_hom/V::
        
                        ---- ----'                             -
                        \    \         q_j            q_hom   /      1
          phi(r_i)  =    |    |   ---------------  +  -----   |dr ---------
                        /    /    |r_i - r_j + l|       V     |   |r - r_i|
                        ---- ----                             / 
                         j    l
        
        r_B : matrix with the lattice basis (in cartesian coordinates).
        
        q_B : probe charges (in units of e).
        
        excludefroml0 : integer specifying if a point charge is not to
        be included in the central (l=0) unit cell. Used for Madelung
        constants.
        """
        E0 = 0.
        if excludefroml0 is None:
            excludefroml0 = np.zeros(len(q_B), dtype=int)
        if excludefroml0 in range(len(q_B)):
            i = excludefroml0
            excludefroml0 = np.zeros(len(q_B), dtype=int)
            excludefroml0[i] = 1
        
        assert sum(excludefroml0) <= 1

        for i, q in enumerate(q_B): # potential from point charges
            rprime = r - r_B[i]
            absr = np.linalg.norm(rprime)
            E0 += q * self.get_sum_real_ij(rprime)
            E0 += q * self.get_sum_recip_ij(rprime)
            if excludefroml0[i]: # if sum over l not 0
                if absr < 1e-14:
                    # lim r -> 0 erf(r G) / r = 2 * G / sqrt(pi)
                    E0 -= q * 2. * self.G / np.sqrt(pi)
                else:
                    E0 -= q * erf(absr * self.G) / absr
            else: # if sum over all l
                E0 += q * erfc(absr * self.G) / absr

        # correct for compensating homogeneous background
        q_hom = -sum(q_B)
        E0 += q_hom * pi / (self.G**2 * self.Vcell)
        
        return E0


def madelung(cell):
    return Ewald(cell).get_electrostatic_potential(np.zeros(3), np.zeros(3),
                                                   [-1], 0)
