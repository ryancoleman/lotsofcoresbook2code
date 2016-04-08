"""Helper classes for doing jellium calculations."""

from math import pi

import numpy as np
from ase import Atoms
from ase.units import Bohr

from gpaw.poisson import PoissonSolver


class JelliumPoissonSolver(PoissonSolver):
    """Jellium Poisson solver."""
    
    mask_g = None  # where to put the jellium
    rs = None  # Wigner Seitz radius
    
    def get_mask(self, r_gv):
        """Choose which grid points are inside the jellium.

        r_gv: 4-dimensional ndarray
            positions of the grid points in Bohr units.

        Return ndarray of ones and zeros indicating where the jellium
        is.  This implementation will put the positive background in the
        whole cell.  Overwrite this method in subclasses."""
        
        return self.gd.zeros() + 1.0
        
    def initialize(self):
        PoissonSolver.initialize(self)
        r_gv = self.gd.get_grid_point_coordinates().transpose((1, 2, 3, 0))
        self.mask_g = self.get_mask(r_gv).astype(float)
        self.volume = self.gd.comm.sum(self.mask_g.sum()) * self.gd.dv
        
    def solve(self, phi, rho, eps=None, charge=0, maxcharge=1e-6,
              zero_initial_phi=False):

        if eps is None:
            eps = self.eps
        
        self.rs = (3 / pi / 4 * self.volume / charge)**(1 / 3.0)
        
        rho -= self.mask_g * (charge / self.volume)
        niter = self.solve_neutral(phi, rho, eps=eps)
        return niter


class JelliumSurfacePoissonSolver(JelliumPoissonSolver):
    def __init__(self, z1, z2, **kwargs):
        """Put the positive background charge where z1 < z < z2.

        z1: float
            Position of lower surface in Angstrom units.
        z2: float
            Position of upper surface in Angstrom units."""
        
        PoissonSolver.__init__(self, **kwargs)
        self.z1 = (z1 - 0.0001)  / Bohr
        self.z2 = (z2 - 0.0001) / Bohr

    def get_mask(self, r_gv):
        return np.logical_and(r_gv[:, :, :, 2] > self.z1,
                              r_gv[:, :, :, 2] < self.z2)
