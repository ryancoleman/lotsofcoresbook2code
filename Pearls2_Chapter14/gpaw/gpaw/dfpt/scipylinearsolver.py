import numpy as np
import scipy.sparse.linalg as sla

class ScipyLinearSolver:
    """Wrapper class for the linear solvers in ``scipy.sparse.linalg``.

    Solve the linear system of equations

            A * x = b

    where A is a linear operator and b is the known rhs. The linear operator
    provided as argument in the ``solve`` method must have a ``shape``
    attribute (a tuble (M,N) where M and N give the size of the corresponding
    matrix) and a ``dtype`` attribute giving datatype of the matrix entries.
    
    """
    
    # Supported solvers
    solvers = {'cg':       sla.cg,       # symmetric positive definite 
               'minres':   sla.minres,   # symmetric indefinite
               'gmres':    sla.gmres,    # non-symmetric
               'bicg':     sla.bicg,     # non-symmetric
               'cgs':      sla.cgs,      # similar to bicg
               'bicgstab': sla.bicgstab, # 
               'qmr':      sla.qmr,      #
               # 'lgmres': sla.lgmres, # scipy v. 0.8.0
               }
    
    def __init__(self, method='gmres', preconditioner=None, tolerance=1e-5,
                 max_iter=100):
        """Initialize the linear solver.

        Parameters
        ----------
        method: str
            One of the supported linear solvers in scipy.
        preconditioner: LinearOperator
            Instance of 1) class ``LinearOperator`` from the
            ``scipy.sparse.linalg`` package, or 2) user defined class with
            required attributes.

        """

        if method not in ScipyLinearSolver.solvers:
            raise RuntimeError("Unsupported solver %s" % method)
                                   
        self.solver = ScipyLinearSolver.solvers[method]
        self.pc = preconditioner

        self.tolerance = tolerance
        self.max_iter = max_iter

        # Iteration counter
        self.i = None
        
    def solve(self, A, x_nG, b_nG):
        """Solve linear system Ax = b.

        Parameters
        ----------
        A: LinearOperator
            Instance of 1) class ``LinearOperator`` from the
            ``scipy.sparse.linalg`` package, or 2) user defined class with
            required attributes.
        x_nG: ndarray
            Vector to be solved for (must contain initial starting vector on
            input).
        b_nG: ndarray
            Vector with the right-hand side of the linear equation.
            
        """

        assert x_nG.shape == b_nG.shape
        shape = x_nG.shape

        # Initialize iteration counter
        self.i = 0
        
        # Reshape arrays for scipy
        x_0 = x_nG.ravel()
        b = b_nG.ravel()
        
        x, info = self.solver(A, b, x0=x_0, maxiter=self.max_iter, M=self.pc,
                              tol=self.tolerance, callback=self.iteration)

        x_nG[:] = x.reshape(shape)
        
        return self.i, info
            
    def iteration(self, x_i):
        """Passed as callback function to the scipy-routine."""

        # Increment iterator counter
        self.i += 1
