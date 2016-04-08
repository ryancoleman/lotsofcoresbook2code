from ase.lattice.surface import fcc111
atoms = fcc111('Al', size=(1, 1, 2))
atoms.center(vacuum=4.0, axis=2)

from gpaw import GPAW
calc = GPAW(mode='pw',
            kpts=(4, 4, 1),
            symmetry='off',
            txt='al111.txt')
atoms.set_calculator(calc)
energy = atoms.get_potential_energy()
calc.write('al111.gpw', 'all')
