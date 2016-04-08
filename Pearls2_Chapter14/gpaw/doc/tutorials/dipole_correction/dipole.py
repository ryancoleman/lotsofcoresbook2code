from ase.lattice.surface import fcc100, add_adsorbate
from gpaw import GPAW
from gpaw.poisson import PoissonSolver
from gpaw.dipole_correction import DipoleCorrection

slab = fcc100('Al', (2, 2, 2), a=4.05, vacuum=7.5)
add_adsorbate(slab, 'Na', 4.0)
slab.center(axis=2)

slab.calc = GPAW(txt='zero.txt',
                 xc='PBE',
                 setups={'Na': '1'},
                 kpts=(4, 4, 1))
e1 = slab.get_potential_energy()
slab.calc.write('zero.gpw')

slab.pbc = True
slab.calc.set(txt='periodic.txt')
e2 = slab.get_potential_energy()
slab.calc.write('periodic.gpw')

slab.pbc = (True, True, False)
slab.calc.set(poissonsolver=DipoleCorrection(PoissonSolver(), 2),
              txt='corrected.txt')
e3 = slab.get_potential_energy()
slab.calc.write('corrected.gpw')
