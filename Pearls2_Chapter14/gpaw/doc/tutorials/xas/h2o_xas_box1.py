import numpy as np
from ase.structure import molecule
from gpaw import GPAW, setup_paths
setup_paths.insert(0, '.')

atoms = molecule('H2O')

h = 0.2

for L in np.arange(4, 14, 2) * 8 * h:
    atoms.set_cell((L, L, L))
    atoms.center()
    calc = GPAW(xc='PBE',
                h=h,
                nbands=-40,
                eigensolver='cg',
                setups={'O': 'hch1s'})
    atoms.set_calculator(calc)
    e1 = atoms.get_potential_energy()
    calc.write('h2o_hch_%.1f.gpw' % L)
