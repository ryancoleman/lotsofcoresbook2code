"""This module implements the linear operator in the Sternheimer equation.

Sternheimer equation::

       /           \             q
       | H - eps   | P      |dPsi  > = - P      dV  |Psi  > ,
       \        nk /  c,k+q      nk       c,k+q   q     nk
   
where P_c is the projection operator onto the unoccupied states. For a
perturbation with a specific q-vector, only the projections onto states at
k+q are needed.

"""

import numpy as np

from gpaw.utilities import unpack
from gpaw.utilities.blas import gemm
from gpaw.fd_operators import Laplace

class SternheimerOperator:
    """Class implementing the linear operator in the Sternheimer equation.
    
    The main purpose of this class is to implement the multiplication of the
    linear operator (H - eps_nk) with a vector in the ``apply`` member
    function. An additional ``matvec`` member function has been defined so that
    instances of this class can be passed as a linear operator to scipy's
    iterative Krylov solvers in ``scipy.sparse.linalg``.

    """
    
    def __init__(self, hamiltonian, wfs, gd, dtype=float):
        """Save useful objects for the Sternheimer operator.

        Parameters
        ----------
        hamiltonian: Hamiltonian
            Hamiltonian for a ground-state calculation.
        wfs: WaveFunctions
            Ground-state wave-functions.
        gd: GridDescriptor
            Grid on which the operator is defined.
        dtype: dtype
            dtype of the wave-function being operated on.

        """

        self.hamiltonian = hamiltonian
        self.kin = Laplace(gd, scale=-0.5, n=3, dtype=dtype)
        self.kpt_u = wfs.kpt_u
        self.pt = wfs.pt
        self.gd = gd

        # Variables for k-point and band index
        self.k = None
        self.n = None
        self.kplusq = None
        
        # For scipy's linear solver
        N = np.prod(gd.n_c)
        self.shape = (N, N)
        self.dtype = dtype

    def set_k(self, k):
        """Set k-vector for the Bloch-state in consideration."""

        self.k = k
        
    def set_band(self, n):
        """Set band index for the Bloch-state in consideration.

        Note that n different from None has implications for the interpretation
        of the input vector in the ``apply`` member function.

        """

        self.n = n

    def set_kplusq(self, kplusq):
        """Set q-vector of the perturbing potential.

        Parameters
        ----------
        kplusq: int
           Index of the k+q vector.

        """

        self.kplusq = kplusq

    def apply(self, x_nG, y_nG):
        """Apply the linear operator in the Sternheimer equation to a vector.

        For the eigenvalue term the k-point is the one of the state.
        For the other terms the k-point to be used is the one given by the k+q
        phase of the first-order of the state. Only for q=0 do the two coincide.
        
        Parameters
        ----------
        x_nG: ndarray
            Vector(s) to which the Sternheimer operator is applied.
        y_nG: ndarray
            Resulting vector(s).
            
        """

        assert x_nG.ndim in (3, 4)
        assert x_nG.shape == y_nG.shape
        assert tuple(self.gd.n_c) == x_nG.shape[-3:]
        assert self.k is not None

        # k
        kpt = self.kpt_u[self.k]
        # k+q
        kplusqpt = self.kpt_u[self.kplusq]
        
        # Kintetic energy
        # k+q
        self.kin.apply(x_nG, y_nG, phase_cd=kplusqpt.phase_cd)

        # Local part of effective potential - no phase !!
        self.hamiltonian.apply_local_potential(x_nG, y_nG, kpt.s)
        
        # Non-local part from projectors (coefficients can not be reused)
        shape = x_nG.shape[:-3]
        P_ani = self.pt.dict(shape=shape)

        # k+q 
        self.pt.integrate(x_nG, P_ani, q=kplusqpt.k)

        for a, P_ni in P_ani.items():
            dH_ii = unpack(self.hamiltonian.dH_asp[a][kpt.s])
            P_ani[a] = np.dot(P_ni, dH_ii)
        # k+q
        self.pt.add(y_nG, P_ani, q=kplusqpt.k)

        # XXX Eigenvalue term
        if self.n is not None:
            # k
            y_nG -= kpt.eps_n[self.n] * x_nG
        else:
            for n, x_G in enumerate(x_nG):
                # k
                y_nG[n] -= kpt.eps_n[n] * x_G

        # Project out undesired (numerical) components
        # k+q
        self.project(y_nG)

    def project(self, x_nG):
        """Project the vector onto the unoccupied states at k+q.

        The projection operator is defined as follows::

                      --                    --             
             P      = >  |Psi ><Psi | = 1 - >  |Psi ><Psi |
              c,k+q   --     c     c        --     v     v 
                       c    k+q   k+q        v    k+q   k+q
        
        """

        # It might be a good idea to move this functionality to its own class
        assert x_nG.ndim == 4
        assert self.kplusq is not None

        # k+q
        kplusqpt = self.kpt_u[self.kplusq]
        # Occupied wave function
        psit_nG = kplusqpt.psit_nG

        # Project out occupied sub-space using blas routine gemm
        m = len(x_nG)
        n = len(psit_nG)
        proj_mn = np.zeros((m, n), dtype=self.dtype)
        gemm(self.gd.dv, psit_nG, x_nG, 0.0, proj_mn, 'c')
        gemm(-1.0, psit_nG, proj_mn, 1.0, x_nG)

        # Project out one orbital at a time
        ## for n, psit_G in enumerate(psit_nG):
        ## 
        ##     proj = self.gd.integrate(psit_G.conjugate() * x_nG)
        ##     x_nG -= proj * psit_G
        
    def matvec(self, x):
        """Matrix-vector multiplication for scipy's Krylov solvers.

        This is a wrapper around the ``apply`` member function above. It allows
        to multiply the sternheimer operator onto the 1-dimensional scipy
        representation of a gpaw grid vector(s).
        
        Parameters
        ----------
        a: ndarray
            1-dimensional array holding the representation of a gpaw grid
            vector (possibly a set of vectors).

        """

        # Check that a is 1-dimensional
        assert len(x.shape) == 1
        
        # Find the number of states in x
        grid_shape = tuple(self.gd.n_c)
        assert ((x.size % np.prod(grid_shape)) == 0), ("Incompatible array " +
                                                       "shapes")
        # Number of states in the vector
        N = x.size / np.prod(grid_shape)
        
        # Output array
        y_nG = self.gd.zeros(n=N, dtype=self.dtype)
        shape = y_nG.shape

        assert x.size == np.prod(y_nG.shape)
        
        x_nG = x.reshape(shape)

        self.apply(x_nG, y_nG)
        
        y = y_nG.ravel()
        
        return y
