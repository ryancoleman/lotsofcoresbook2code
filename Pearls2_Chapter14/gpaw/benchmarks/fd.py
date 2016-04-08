""" Test the finite difference stencil """

from gpaw.grid_descriptor import GridDescriptor
from gpaw.fd_operators import Laplace
from gpaw.test import equal
import numpy as np

gpts = 64
nbands = 120
gd = GridDescriptor([gpts, gpts, gpts])

a = gd.empty(nbands)
b = gd.empty(nbands)

np.random.seed(10)
a[:] = np.random.random(a.shape)

op = Laplace(gd, 1.0, 3).apply
op(a, b)

equal(np.sum(b), -7.43198208966e-05, 1e-9)
