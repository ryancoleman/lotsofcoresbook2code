import numpy as np
from ase import Atoms
from gpaw import GPAW
from gpaw.mixer import FFTMixer
from gpaw.wavefunctions.pw import PW
from gpaw.test import equal

bulk = Atoms('Li', pbc=True)
k = 4
calc = GPAW(mode=PW(200), kpts=(k, k, k), mixer=FFTMixer(), eigensolver='rmm-diis')
bulk.set_calculator(calc)
bulk.set_cell((2.6, 2.6, 2.6))
e = bulk.get_potential_energy()
nit = calc.get_number_of_iterations()
equal(e, -1.98481281259, 1.0e-8)
