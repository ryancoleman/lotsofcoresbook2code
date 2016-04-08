"""Wigner-Seitz truncated coulomb interaction.

See:

    Ravishankar Sundararaman and T. A. Arias:
    Phys. Rev. B 87, 165122 (2013)
    
    Regularization of the Coulomb singularity in exact exchange by
    Wigner-Seitz truncated interactions: Towards chemical accuracy
    in nontrivial systems
"""

import sys
from math import pi

import numpy as np
from ase.units import Bohr
from ase.utils import prnt

import gpaw.mpi as mpi
from gpaw.utilities import erf
from gpaw.fftw import get_efficient_fft_size
from gpaw.grid_descriptor import GridDescriptor


class WignerSeitzTruncatedCoulomb:
    def __init__(self, cell_cv, nk_c, txt=sys.stdout):
        self.nk_c = nk_c
        bigcell_cv = cell_cv * nk_c[:, np.newaxis]
        L_c = (np.linalg.inv(bigcell_cv)**2).sum(0)**-0.5
        
        rc = 0.5 * L_c.min()
        prnt('Inner radius for %dx%dx%d Wigner-Seitz cell: %.3f Ang' %
             (tuple(nk_c) + (rc * Bohr,)), file=txt)
        
        self.a = 5 / rc
        prnt('Range-separation parameter: %.3f Ang^-1' % (self.a / Bohr),
             file=txt)
        
#        nr_c = [get_efficient_fft_size(2 * int(L * self.a * 1.5))
        nr_c = [get_efficient_fft_size(2 * int(L * self.a * 3.0))
                for L in L_c]
        prnt('FFT size for calculating truncated Coulomb: %dx%dx%d' %
             tuple(nr_c), file=txt)
        
        self.gd = GridDescriptor(nr_c, bigcell_cv, comm=mpi.serial_comm)
        v_R = self.gd.empty()
        v_i = v_R.ravel()
        
        pos_iv = self.gd.get_grid_point_coordinates().reshape((3, -1)).T
        corner_jv = np.dot(np.indices((2, 2, 2)).reshape((3, 8)).T, bigcell_cv)
        for i, pos_v in enumerate(pos_iv):
            r = ((pos_v - corner_jv)**2).sum(axis=1).min()**0.5
            if r == 0:
                v_i[i] = 2 * self.a / pi**0.5
            else:
                v_i[i] = erf(self.a * r) / r
                
        self.K_Q = np.fft.fftn(v_R) * self.gd.dv
        
    def get_potential(self, pd):
        q_c = pd.kd.bzk_kc[0]
        shift_c = (q_c * self.nk_c).round().astype(int)
        max_c = self.gd.N_c // 2
        K_G = pd.zeros()
        N_c = pd.gd.N_c
        for G, Q in enumerate(pd.Q_qG[0]):
            Q_c = (np.unravel_index(Q, N_c) + N_c // 2) % N_c - N_c // 2
            Q_c = Q_c * self.nk_c + shift_c
            if (abs(Q_c) < max_c).all():
                K_G[G] = self.K_Q[tuple(Q_c)]

        G2_G = pd.G2_qG[0]
        a = self.a
        if pd.kd.gamma:
            K_G[0] += pi / a**2
        else:
            K_G[0] += 4 * pi * (1 - np.exp(-G2_G[0] / (4 * a**2))) / G2_G[0]
        K_G[1:] += 4 * pi * (1 - np.exp(-G2_G[1:] / (4 * a**2))) / G2_G[1:]
        assert pd.dtype == complex
        return K_G
