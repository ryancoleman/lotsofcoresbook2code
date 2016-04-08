# {   -1,   -1}:  On entry to PDLASRT parameter number    9 had an illegal value

# works with 'sl_default': (2, 2, 32)

from ase.lattice.surface import fcc100, add_adsorbate
from gpaw import GPAW, ConvergenceError
from gpaw.mpi import world
from gpaw.utilities import compiled_with_sl

assert world.size == 4

slab = fcc100('Cu', size=(2, 2, 2))
add_adsorbate(slab,'O', 1.1, 'hollow')
slab.center(vacuum=3.0, axis=2)

calc = GPAW(mode='lcao',
            kpts=(2, 2, 1),
            txt='-',
            maxiter=1,)

if compiled_with_sl():
    calc.set(parallel={'domain': (1, 1, 4), 'sl_default': (2, 2, 64)})

slab.set_calculator(calc)
try:
    slab.get_potential_energy()
except ConvergenceError:
    pass
