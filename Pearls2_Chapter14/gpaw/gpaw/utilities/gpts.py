import numpy as np
from ase.units import Bohr, Hartree

from gpaw.utilities import h2gpts
from gpaw.wavefunctions.pw import PW
from gpaw.fftw import get_efficient_fft_size
from gpaw.wavefunctions.fd import FD


def get_number_of_grid_points(cell_cv, h=None, mode=None, realspace=None,
                              symmetry=None):
    if mode is None:
        mode = FD()

    if realspace is None:
        realspace = mode.name != 'pw'

    if h is None:
        if mode.name == 'pw':
            h = np.pi / (4 * mode.ecut)**0.5
        elif mode.name == 'lcao' and not realspace:
            h = np.pi / (4 * 340 / Hartree)**0.5
        else:
            h = 0.2 / Bohr

    if realspace or mode.name == 'fd':
        N_c = h2gpts(h, cell_cv, 4)
    else:
        N_c = h2gpts(h, cell_cv, 1)
        if symmetry is None:
            N_c = np.array([get_efficient_fft_size(N) for N in N_c])
        else:
            N_c = np.array([get_efficient_fft_size(N, n)
                            for N, n in zip(N_c, symmetry.gcd_c)])
    
    return N_c
