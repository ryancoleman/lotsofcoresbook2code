# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

import numpy as np

from gpaw import debug
from gpaw.utilities import is_contiguous
import _gpaw


class Spline:
    def __init__(self, l, rmax, f_g):
        """Spline(l, rcut, list) -> Spline object

        The integer l gives the angular momentum quantum number and
        the list contains the spline values from r=0 to r=rcut.

        The array f_g gives the radial part of the function on the grid.
        The radial function is multiplied by a real solid spherical harmonics
        (r^l * Y_lm).
        """
        assert 0.0 < rmax
        f_g = np.array(f_g, float)
        # Copy so we don't change the values of the input array
        f_g[-1] = 0.0
        self.spline = _gpaw.Spline(l, rmax, f_g)

    def get_cutoff(self):
        """Return the radial cutoff."""
        return self.spline.get_cutoff()

    def get_angular_momentum_number(self):
        """Return the angular momentum quantum number."""
        return self.spline.get_angular_momentum_number()

    def get_value_and_derivative(self, r):
        """Return the value and derivative."""
        return self.spline.get_value_and_derivative(r)

    def __call__(self, r):
        assert r >= 0.0
        return self.spline(r)

    def map(self, r_x):
        return np.vectorize(self, [float])(r_x)

    def get_functions(self, gd, start_c, end_c, spos_c):
        h_cv = gd.h_cv
        # start_c is the new origin so we translate gd.beg_c to start_c
        origin_c = np.array([0,0,0])
        pos_v = np.dot(spos_c, gd.cell_cv) - np.dot(start_c, h_cv)
        A_gm, G_b = _gpaw.spline_to_grid(self.spline, origin_c, end_c-start_c,
                                         pos_v, h_cv, end_c-start_c, origin_c)

        if debug:
            assert G_b.ndim == 1 and G_b.shape[0] % 2 == 0
            assert is_contiguous(G_b, np.intc)
            assert A_gm.shape[:-1] == np.sum(G_b[1::2]-G_b[::2])

        indices_gm, ng, nm = self.spline.get_indices_from_zranges(start_c,
                                                                  end_c, G_b)
        shape = (nm,) + tuple(end_c-start_c)
        work_mB = np.zeros(shape, dtype=A_gm.dtype)
        np.put(work_mB, indices_gm, A_gm)
        return work_mB


## class rspline:
##     def __init__(self, r_g, f_g, l=0):
##         self.rcut = r_g[-1]
##         self.l = l
##         ...
        
##     def __call__(self, r):
##         return self.spline(r)*r**l
