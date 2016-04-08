from ase import Atoms
from gpaw import GPAW
from gpaw.xc.noncollinear import NonCollinearFunctional, \
     NonCollinearLCAOEigensolver, NonCollinearMixer
from gpaw.xc import XC

h = Atoms('H', magmoms=[1])
h.center(vacuum=2)
h.calc = GPAW(txt='c.txt', mode='lcao', basis='dz(dzp)', h=0.25)
e0 = h.get_potential_energy()

h.set_initial_magnetic_moments()
h.set_initial_magnetic_moments([(0.1, 0.2, 0.3)])
h.calc = GPAW(txt='nc.txt', mode='lcao', basis='dz(dzp)', h=0.25,
              xc=NonCollinearFunctional(XC('LDA')),
              mixer=NonCollinearMixer(),
              eigensolver=NonCollinearLCAOEigensolver())
e = h.get_potential_energy()

assert abs(e - e0) < 2e-5, (e, e0)
from ase.io import write
write('h.traj', h)

