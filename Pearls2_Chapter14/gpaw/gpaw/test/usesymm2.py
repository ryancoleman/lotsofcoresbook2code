from math import sqrt
from ase import Atoms
from gpaw import GPAW

a = 2.4437
# Set up graphene structure
vacuum = 3
atoms = Atoms(symbols='C2',
              positions=[(0.5 * a, -sqrt(3) / 6 * a, vacuum),
                         (0.5 * a, +sqrt(3) / 6 * a, vacuum)],
              cell=[(0.5 * a, -0.5 * 3**0.5 * a, 0),
                    (0.5 * a, +0.5 * 3**0.5 * a, 0),
                    (0.0, 0.0, 2 * vacuum)])
atoms.pbc = True

# Gamma-point calculation of graphene
calc = GPAW(h=0.2, width=0.15, kpts=(1, 1, 1), xc='LDA')
atoms.set_calculator(calc)
atoms.get_potential_energy()


kpts = [(1 / 2.0, 1 / 3.0, 0)]

# Calculate one K-point with symmetry:
calc.set(kpts=kpts, fixdensity=True)
calc.get_potential_energy()
eigs_True = calc.get_eigenvalues(kpt=0)

# Redo with the same K-point without symmetry:
calc.set(kpts=kpts,
         symmetry='off',
         fixdensity=True)
calc.get_potential_energy()
eigs_False = calc.get_eigenvalues(kpt=0)

assert abs(eigs_True[0] - eigs_False[0]) < 1e-4
