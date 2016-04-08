import numpy as np

from ase.units import Bohr
from ase.parallel import paropen
from ase.utils import prnt

from gpaw.spherical_harmonics import Y
from gpaw.utilities.tools import coordinates


class Multipole:

    """Expand a function on the grid in multipole moments
    relative to a given center.

    center: Vector [Angstrom]
    """

    def __init__(self, center, calculator=None, lmax=6):

        self.center = center / Bohr
        self.lmax = lmax

        self.gd = None
        self.y_Lg = None
        self.l_L = None

        if calculator is not None:
            self.initialize(calculator.density.finegd)

    def initialize(self, gd):
        """Initialize Y_L arrays"""

        self.gd = gd

        r_cg, r2_g = coordinates(gd, self.center, tiny=1.e-78)
        r_g = np.sqrt(r2_g)
        rhat_cg = r_cg / r_g

        self.l_L = []
        self.y_Lg = []
        npY = np.vectorize(Y, (float,), 'spherical harmonic')
        L = 0
        for l in range(self.lmax + 1):
            for m in range(2 * l + 1):
                self.y_Lg.append(
                    np.sqrt(4 * np.pi / (2 * l + 1)) * r_g ** l *
                    npY(L, rhat_cg[0], rhat_cg[1], rhat_cg[2])
                )
                self.l_L.append(l)
                L += 1

    def expand(self, f_g):
        """Expand a function f_g in multipole moments
        units [e * Angstrom**l]"""

        assert(f_g.shape == self.gd.empty().shape)

        q_L = []
        for L, y_g in enumerate(self.y_Lg):
            q_L.append(self.gd.integrate(f_g * y_g))
            q_L[L] *= Bohr ** self.l_L[L]

        return np.array(q_L)

    def to_file(self, calculator,
                filename='multipole.dat',
                mode='a'):
        """Expand the charge distribution in multipoles and write
        the result to a file"""

        if self.gd is None:
            self.initialize(calculator.density.finegd)
        q_L = self.expand(-calculator.density.rhot_g)

        f = paropen(filename, mode)

        prnt('# Multipole expansion of the charge density', file=f)
        prnt('# center =', self.center * Bohr, 'Angstrom', file=f)
        prnt('# lmax =', self.lmax, file=f)
        prnt(('# see https://trac.fysik.dtu.dk/projects/gpaw/browser/' +
              'trunk/c/bmgs/sharmonic.py'), file=f)
        prnt('# for the definition of spherical harmonics', file=f)
        prnt('# l  m    q_lm[|e| Angstrom**l]', file=f)

        L = 0
        for l in range(self.lmax + 1):
            for m in range(-l, l + 1):
                prnt('{0:2d} {1:3d} {2:g}'.format(l, m, q_L[L]), file=f)
                L += 1
        f.close()
