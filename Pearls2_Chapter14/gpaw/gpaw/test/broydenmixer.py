from ase import Atom, Atoms
from gpaw import GPAW
from gpaw.mixer import BroydenMixer
from gpaw.test import equal

a = 2.7
bulk = Atoms([Atom('Li')], pbc=True, cell=(a, a, a))
k = 2
g = 16
calc = GPAW(gpts=(g, g, g), kpts=(k, k, k), nbands=2,
                  mixer=BroydenMixer())
bulk.set_calculator(calc)
e = bulk.get_potential_energy()
niter = calc.get_number_of_iterations()
calc.write('Li.gpw')
calc2 = GPAW('Li.gpw')

energy_tolerance = 0.00005
niter_tolerance = 0
equal(e, -1.20258, energy_tolerance)
