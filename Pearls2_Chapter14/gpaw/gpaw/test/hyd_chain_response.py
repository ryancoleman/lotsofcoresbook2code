# Test of the handling of degenerate bands in response code
from ase import Atoms
from gpaw import GPAW, PW
from gpaw.response.df import DielectricFunction
from gpaw.test import findpeak, equal
import numpy as np


def get_hydrogen_chain_dielectric_function(NH, NK):
    a = Atoms('H', cell=[1, 1, 1], pbc=True)
    a.center()
    a = a.repeat((1, 1, NH))
    a.calc = GPAW(mode=PW(200), kpts={'size': (1, 1, NK), 'gamma': True},
                  parallel={'band': 1}, dtype=complex, gpts=(10, 10, 10 * NH))
    a.get_potential_energy()
    a.calc.diagonalize_full_hamiltonian(nbands=2 * NH)
    a.calc.write('H_chain.gpw', 'all')

    DF = DielectricFunction('H_chain.gpw', ecut=1e-3, hilbert=False,
                            omega2=np.inf, intraband=False)
    eps_NLF, eps_LF = DF.get_dielectric_function(direction='z')
    omega_w = DF.get_frequencies()
    return omega_w, eps_LF

NH_i = [2**n for n in [0, 4]]
NK_i = [2**n for n in [6, 2]]

opeak_old = np.nan
peak_old = np.nan

for NH, NK in zip(NH_i, NK_i):
    omega_w, eps_w = get_hydrogen_chain_dielectric_function(NH, NK)
    eels_w = -(1. / eps_w).imag
    opeak, peak = findpeak(omega_w, eels_w)

    # Test for consistency
    if not np.isnan(opeak_old):
        equal(opeak, opeak_old, tolerance=1e-3)
        equal(peak, peak_old, tolerance=1e-3)
    opeak_old = opeak
    peak_old = peak
