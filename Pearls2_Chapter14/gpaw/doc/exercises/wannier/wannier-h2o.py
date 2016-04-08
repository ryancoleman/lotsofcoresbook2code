from ase.structure import molecule
from ase.visualize import view
from ase import Atoms
from gpaw import GPAW
from gpaw.wannier import Wannier

calc = GPAW(nbands=4)

atoms = molecule('H2O')
atoms.center(vacuum=3.)
atoms.set_calculator(calc)
atoms.get_potential_energy()

# Initialize the Wannier class
w = Wannier(calc)
w.localize()
centers = w.get_centers()

view(atoms + Atoms(symbols='X4', positions=centers))
