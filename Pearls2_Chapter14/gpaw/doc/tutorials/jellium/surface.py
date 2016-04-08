import numpy as np
from ase import Atoms
from ase.units import Bohr
from gpaw.jellium import JelliumSurfacePoissonSolver
from gpaw import GPAW, Mixer

rs = 5.0 * Bohr  # Wigner-Seitz radius
h = 0.2          # grid-spacing
a = 8 * h        # lattice constant
v = 3 * a        # vacuum
L = 10 * a       # thickness
k = 12           # number of k-points (k*k*1)

ne = a**2 * L / (4 * np.pi / 3 * rs**3)

ps = JelliumSurfacePoissonSolver(z1=v, z2=v + L)
surf = Atoms(pbc=(True, True, False),
             cell=(a, a, v + L + v))
surf.calc = GPAW(poissonsolver=ps,
                 xc='LDA_X+LDA_C_WIGNER',
                 eigensolver='dav',
                 charge=-ne,
                 kpts=[k, k, 1],
                 h=h,
                 maxiter=300,
                 convergence={'density': 0.001},
                 mixer=Mixer(0.03, 7, 100),
                 nbands=int(ne / 2) + 15,
                 txt='surface.txt')
e = surf.get_potential_energy()
surf.calc.write('surface.gpw')
