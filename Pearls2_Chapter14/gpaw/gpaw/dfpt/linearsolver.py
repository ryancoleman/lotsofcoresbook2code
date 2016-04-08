import numpy as np

from scipy.sparse.linalg import cgs

class LinearSolver:
    
    def __init__(self, method=None, preconditioner=None, tolerance = 1e-10,
                 max_iter = 10000):
        """Init method."""

        # self.solver = solvers[method]
        # self.pc = pcs[preconditioner]
        self.tolerance = tolerance
        self.max_iter = max_iter
        
    def solve(self, A, x_G, b_G):
        """Solve linear system A * x = b."""

       
        ##########################################
        # Simple steepest descent implementation #
        ##########################################

        # Residue: r = A*x - b
        r_G = np.zeros_like(x_G)

        A.dot(x_G,r_G)
        # A.project(r_G)
        r_G -= b_G

        for iter in range(self.max_iter):
         
            x_G -= 0.01 * r_G
            A.dot(x_G,r_G)
            # Try and project out the undersired components in the residue
            # A.project(r_G)
            r_G -= b_G
            error = A.norm(r_G)

            if error < self.tolerance:
                break

        return iter
            
    def apply_preconditioner(self, x, b):
        """Solves preconditioner equation.

        Parameters
        ----------
        psi: List of coarse grids
            the known wavefunctions
        psin: List of coarse grids
            the result

        """

        #self.timer.start('Solve Sternheimer preconditioner')
        if self.preconditioner is not None:
            #self.preconditioner.apply(self.kpt, psi, psin)
            pass
        else:
            psin[:] = psi
        #self.timer.stop('Solve Sternheimer preconditioner')
