from ase import Atoms
from gpaw import GPAW
from gpaw.test import equal


# Self-consistent calculation:
a = 2.5
slab = Atoms('Li', cell=(a, a, 2 * a), pbc=1)
slab.calc = GPAW(kpts=(3,3,1), txt='li.txt',
                 parallel=dict(kpt=1))
slab.get_potential_energy()
slab.calc.write('Li.gpw')

# Gamma point:
e1 = slab.calc.get_eigenvalues(kpt=0)[0]

# Fix density and continue:
kpts = [(0,0,0)]
slab.calc.set(fixdensity=True,
              nbands=5,
              kpts=kpts,
              symmetry='off',
              eigensolver='cg')
slab.get_potential_energy()
e2 = slab.calc.get_eigenvalues(kpt=0)[0]

# Start from gpw-file:
calc = GPAW('Li.gpw', 
            fixdensity=True,
            nbands=5,
            kpts=kpts,
            symmetry='off',
            eigensolver='cg')
calc.scf.reset()
calc.get_potential_energy()
e3 = calc.get_eigenvalues(kpt=0)[0]

equal(e1, e2, 3e-5)
equal(e1, e3, 3e-5)
