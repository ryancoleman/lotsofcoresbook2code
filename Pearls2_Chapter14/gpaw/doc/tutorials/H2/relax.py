from __future__ import print_function

from gpaw import restart
from ase.parallel import paropen as open
from ase.optimize import QuasiNewton


molecule, calc = restart('H2.gpw', txt='H2-relaxed.txt')

e2 = molecule.get_potential_energy()
d0 = molecule.get_distance(0, 1)

fd = open('optimization.txt', 'w')
print('experimental bond length:', file=fd)
print('hydrogen molecule energy: %5.2f eV' % e2, file=fd)
print('bondlength              : %5.2f Ang' % d0, file=fd)

# Find the theoretical bond length:
relax = QuasiNewton(molecule, logfile='qn.log')
relax.run(fmax=0.05)

e2 = molecule.get_potential_energy()
d0 = molecule.get_distance(0, 1)

print(file=fd)
print('PBE energy minimum:', file=fd)
print('hydrogen molecule energy: %5.2f eV' % e2, file=fd)
print('bondlength              : %5.2f Ang' % d0, file=fd)
fd.close()
