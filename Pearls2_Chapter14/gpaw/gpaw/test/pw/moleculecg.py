# fails with On entry to ZGEMV parameter number 8 had an illegal value

from ase.structure import molecule
from gpaw import GPAW
from gpaw.wavefunctions.pw import PW

m = molecule('H')
m.center(vacuum=2.0)
m.set_calculator(GPAW(mode=PW(), eigensolver='cg'))
m.get_potential_energy()
