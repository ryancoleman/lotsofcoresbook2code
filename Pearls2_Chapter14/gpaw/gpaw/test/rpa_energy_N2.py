from __future__ import print_function
from ase import *
from ase.structure import molecule
from gpaw import *
from gpaw.test import equal
from gpaw.xc.rpa import RPACorrelation
from gpaw.xc.exx import EXX
import numpy as np

ecut = 25

N2 = molecule('N2')
N2.center(vacuum=2.0)

calc = GPAW(mode='pw', dtype=complex, xc='PBE', eigensolver='rmm-diis')
N2.set_calculator(calc)
E_n2_pbe = N2.get_potential_energy()

calc.diagonalize_full_hamiltonian(nbands=104, scalapack=True)
calc.write('N2.gpw', mode='all')

exx = EXX('N2.gpw')
exx.calculate()
E_n2_hf = exx.get_total_energy()

rpa = RPACorrelation('N2.gpw', nfrequencies=8)
E_n2_rpa = rpa.calculate(ecut=[ecut])
                                    
#-------------------------------------------------------------------------

N = molecule('N')
N.set_cell(N2.cell)

calc = GPAW(mode='pw', dtype=complex, xc='PBE', eigensolver='rmm-diis')
N.set_calculator(calc)
E_n_pbe = N.get_potential_energy()

calc.diagonalize_full_hamiltonian(nbands=104, scalapack=True)
calc.write('N.gpw', mode='all')

exx = EXX('N.gpw')
exx.calculate()
E_n_hf = exx.get_total_energy()

rpa = RPACorrelation('N.gpw', nfrequencies=8)
E_n_rpa = rpa.calculate(ecut=[ecut])

print('Atomization energies:')
print('PBE: ', E_n2_pbe - 2 * E_n_pbe)
print('HF: ',  E_n2_hf - 2 * E_n_hf)
print('HF+RPA: ', E_n2_hf - 2 * E_n_hf + E_n2_rpa[0] - 2 * E_n_rpa[0])

equal(E_n2_rpa - 2*E_n_rpa, -1.72, 0.02)
