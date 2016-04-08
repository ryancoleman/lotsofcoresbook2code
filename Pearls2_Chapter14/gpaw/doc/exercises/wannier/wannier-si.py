from ase import Atoms
from ase.visualize import view
from gpaw import GPAW
from gpaw.wannier import Wannier

calc = GPAW('si.gpw')
atoms = calc.get_atoms()

w = Wannier(calc)
w.localize()
centers = w.get_centers()

watoms = atoms + Atoms(symbols='X16', positions=centers)
view(watoms * (2, 2, 2))
