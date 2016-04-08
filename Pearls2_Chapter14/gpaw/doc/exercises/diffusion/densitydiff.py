import numpy as np
from gpaw import restart

slab, calc = restart('ontop.gpw', txt=None)
AuAl_density = calc.get_pseudo_density()

# Remove gold atom and do a clean slab calculation:
del slab[-1]
slab.get_potential_energy()
Al_density = calc.get_pseudo_density()

# Remove Al atoms and do a calculation for Au only:
slab, calc = restart('ontop.gpw', txt=None)
del slab[:-1]
calc.set(kpts=None)
slab.get_potential_energy()
Au_density = calc.get_pseudo_density()

diff = AuAl_density - Au_density - Al_density
np.save('densitydiff.npy', diff)
