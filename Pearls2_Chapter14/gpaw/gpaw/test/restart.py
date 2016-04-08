from __future__ import print_function
import os
from gpaw import GPAW, restart
from ase import Atoms
from gpaw.test import equal
from math import sqrt
import numpy as np

modes = ['gpw']
try:
    import _gpaw_hdf5
    modes.append('hdf5')
except ImportError:
    pass

d = 3.0
atoms = Atoms('Na3', positions=[( 0, 0, 0),
                              ( 0, 0, d),
                              ( 0, d*sqrt(3./4.), d/2.)],
                   magmoms=[1.0, 1.0, 1.0],
                   cell=(3.5, 3.5, 4.+2/3.),
                   pbc=True)

# Only a short, non-converged calcuation
conv = {'eigenstates': 1.24, 'energy':2e-1, 'density':1e-1}
calc = GPAW(h=0.30, nbands=3,
            setups={'Na': '1'},
            convergence=conv)
atoms.set_calculator(calc)
e0 = atoms.get_potential_energy()
niter0 = calc.get_number_of_iterations()
f0 = atoms.get_forces()
m0 = atoms.get_magnetic_moments()
eig00 = calc.get_eigenvalues(spin=0)
eig01 = calc.get_eigenvalues(spin=1)
# Write the restart file(s)
for mode in modes:
    calc.write('tmp.%s' % mode)

del atoms, calc
# Try restarting from all the files
for mode in modes:
    atoms, calc = restart('tmp.%s' % mode)
    e1 = atoms.get_potential_energy()
    try: # number of iterations needed in restart
        niter1 = calc.get_number_of_iterations()
    except: pass
    f1 = atoms.get_forces()
    m1 = atoms.get_magnetic_moments()
    eig10 = calc.get_eigenvalues(spin=0)
    eig11 = calc.get_eigenvalues(spin=1)
    print(e0, e1)
    equal(e0, e1, 1e-10)
    print(f0, f1)
    for ff0, ff1 in zip(f0, f1):
        err = np.linalg.norm(ff0-ff1)
        assert err <= 1e-10
    print(m0, m1)
    for mm0, mm1 in zip(m0, m1):
        equal(mm0, mm1, 1e-10)
    print('A',eig00, eig10)
    for eig0, eig1 in zip(eig00, eig10):
        equal(eig0, eig1, 1e-10)
    print('B',eig01, eig11)
    for eig0, eig1 in zip(eig01, eig11):
        equal(eig0, eig1, 1e-10)

    # Check that after restart everything is writable
    calc.write('tmp2.%s' % mode)
