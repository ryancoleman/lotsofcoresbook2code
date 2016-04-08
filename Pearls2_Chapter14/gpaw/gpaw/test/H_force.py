from __future__ import print_function
from ase import Atom, Atoms
from ase.calculators.test import numeric_force
from gpaw import GPAW, Mixer, FermiDirac
from gpaw.test import equal

a = 4.0
n = 16
atoms = Atoms([Atom('H', [1.234, 2.345, 3.456])],
                    cell=(a, a, a), pbc=True)
calc = GPAW(nbands=1,
            gpts=(n, n, n),
            txt=None,
            mixer=Mixer(0.25, 3, 1),
            convergence={'energy': 1e-7},
            occupations=FermiDirac(0.0))
atoms.set_calculator(calc)
e1 = atoms.get_potential_energy()
niter1 = calc.get_number_of_iterations()
f1 = atoms.get_forces()[0]
for i in range(3):
    f2i = numeric_force(atoms, 0, i)
    print(f1[i]-f2i)
    equal(f1[i], f2i, 0.00025)

energy_tolerance = 0.00006
force_tolerance = 0.0001
niter_tolerance = 0
equal(e1, -0.531042, energy_tolerance)
f1_ref = [-0.291893, -0.305174, -0.35329]
for i in range(3):
    equal(f1[i], f1_ref[i], force_tolerance)
