from __future__ import print_function

import numpy as np
from ase.utils.bee import get_ensemble_energies
from ase.parallel import paropen as open
from gpaw import GPAW 

atom = GPAW('H.gpw', txt=None).get_atoms()
molecule = GPAW('H2.gpw', txt=None).get_atoms()
e1 = atom.get_potential_energy()
e2 = molecule.get_potential_energy()
ea = 2 * e1 - e2

fd = open('ensemble_energies.txt', 'w')
print('PBE:', ea, 'eV', file=fd)

e1i = get_ensemble_energies(atom)
e2i = get_ensemble_energies(molecule)
eai = 2 * e1i - e2i

n = len(eai)
ea0 = np.sum(eai) / n
sigma = (np.sum((eai - ea0)**2) / n)**0.5
print('Best fit:', ea0, '+-', sigma, 'eV', file=fd)
fd.close()

fd = open('ensemble.dat', 'w')
for e in eai:
    print(e, file=fd)
fd.close()
