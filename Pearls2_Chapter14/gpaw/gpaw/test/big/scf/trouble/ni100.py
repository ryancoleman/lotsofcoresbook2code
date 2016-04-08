ncpus = 4
from ase.lattice.surface import fcc100
from gpaw import GPAW
atoms = fcc100(symbol='Ni', size=(1, 1, 9), a=3.52, vacuum=5.5)
atoms.set_initial_magnetic_moments([0.6] * len(atoms))
atoms.calc = GPAW(h=0.18, kpts=(8, 8, 1), xc='PBE')
