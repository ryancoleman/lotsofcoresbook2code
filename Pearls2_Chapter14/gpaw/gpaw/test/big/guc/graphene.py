import numpy as np
from math import sqrt
from ase import Atoms
from ase.lattice.surface import hcp0001
from gpaw import GPAW
from gpaw.test import equal

# Vacuum and hcp lattice parameter for graphene
d = 4.0
a = 2.4437

calc = GPAW(h=0.15, width=0.1, xc='LDA', nbands=-4, txt='-',
            basis='dzp', convergence={'energy': 1e-5, 'density': 1e-5})

# Calculate potential energy per atom for orthogonal unitcell
atoms = hcp0001('C', a=a/sqrt(3), vacuum=d, size=(3,2,1), orthogonal=True)
del atoms[[1,-1]]
atoms.center(axis=0)
atoms.set_calculator(calc)

kpts_c = np.ceil(50 / np.sum(atoms.get_cell()**2, axis=1)**0.5).astype(int)
kpts_c[~atoms.get_pbc()] = 1
calc.set(kpts=kpts_c)
eppa1 = atoms.get_potential_energy() / len(atoms)
F1_av = atoms.get_forces()
equal(np.abs(F1_av).max(), 0, 5e-3)

# Redo calculation with non-orthogonal unitcell
atoms = Atoms(symbols='C2', pbc=(True,True,False),
              positions=[(a/2,-sqrt(3)/6*a,d), (a/2,sqrt(3)/6*a,d)],
              cell=[(a/2,-sqrt(3)/2*a,0), (a/2,sqrt(3)/2*a,0), (0,0,2*d)])
atoms.set_calculator(calc)

kpts_c = np.ceil(50 / np.sum(atoms.get_cell()**2, axis=1)**0.5).astype(int)
kpts_c[~atoms.get_pbc()] = 1
calc.set(kpts=kpts_c)
eppa2 = atoms.get_potential_energy() / len(atoms)
F2_av = atoms.get_forces()
equal(np.abs(F2_av).max(), 0, 5e-3)

equal(eppa1, eppa2, 1e-3)
