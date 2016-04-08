"""This module wraps the gpaw preconditioner for use with scipy solvers."""

import numpy as np

from gpaw.fd_operators import Laplace
from gpaw.preconditioner import Preconditioner


class ScipyPreconditioner:

    def __init__(self, gd, project, dtype=float):
        """Init the gpaw preconditioner.

        Parameters
        ----------
        gd: GridDescriptor
            Coarse grid
            
        """

        self.project = project
        self.gd = gd
        # K-point for the preconditioner
        self.kpt = None
        
        kin = Laplace(gd, scale=-0.5, n=3, dtype=dtype)
        self.pc = Preconditioner(gd, kin, dtype=dtype)
        
        # For scipy's linear solver
        N = np.prod(gd.n_c)
        self.shape = (N,N)
        self.dtype = dtype

    def set_kpt(self, kpt):
        """Set k-point for ``Transformer`` objects inside the preconditioner.

        Parameters
        ----------
        kpt: KPoint or KPointContainer
            Only requirement is that it must have an ``phase_cd`` attribute.

        """

        self.kpt = kpt
        
    def matvec(self, x):
        """Matrix vector multiplication for ``scipy.sparse.linalg`` solvers.

        Parameters
        ----------
        x: ndarry
            1-dimensional array holding the grid representation of the vector

        """

        # Output array
        y_G = self.gd.zeros(dtype=self.dtype)
        shape = y_G.shape

        size = x.size
        assert size == np.prod(shape)
        
        x_G = x.reshape(shape)

        # Call gpaw preconditioner
        y_G = self.pc(x_G, kpt=self.kpt)

        # Project out undesired (numerical) components
        self.project(y_G)
        
        y = y_G.ravel()
        
        return y
