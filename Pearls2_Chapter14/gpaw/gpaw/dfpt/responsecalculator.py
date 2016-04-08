from __future__ import print_function
"""This module implements a linear response calculator class."""

__all__ = ["ResponseCalculator"]

import numpy as np

from gpaw.transformers import Transformer
from gpaw.dfpt.mixer import BaseMixer
from gpaw.dfpt.sternheimeroperator import SternheimerOperator
from gpaw.dfpt.scipylinearsolver import ScipyLinearSolver
from gpaw.dfpt.preconditioner import ScipyPreconditioner


class ResponseCalculator:
    """This class is a calculator for the sc density variation.

    From the given perturbation, the set of coupled equations for the
    first-order density response is solved self-consistently.

    Parameters
    ----------
    max_iter: int
        Maximum number of iterations in the self-consistent evaluation of
        the density variation.
    tolerance_sc: float
        Tolerance for the self-consistent loop measured in terms of
        integrated absolute change of the density derivative between two
        iterations.
    tolerance_sternheimer: float
        Tolerance for the solution of the Sternheimer equation -- passed to
        the ``LinearSolver``.
    beta: float (0 < beta < 1)
        Mixing coefficient.
    nmaxold: int
        Length of history for the mixer.
    weight: int
        Weight for the mixer metric (=1 -> no metric used).
        
    """

    parameters = {'verbose':               False,
                  'max_iter':              100,
                  'max_iter_krylov':       1000,
                  'krylov_solver':         'cg',
                  'tolerance_sc':          1.0e-5,
                  'tolerance_sternheimer': 1.0e-4,
                  'use_pc':                True,
                  'beta':                  0.1,
                  'nmaxold':               6,
                  'weight':                50
                  }
    
    def __init__(self, calc, wfs, poisson_solver=None, dtype=float, **kwargs):
        """Store calculator etc.

        Parameters
        ----------
        calc: Calculator
            Calculator instance containing a ground-state calculation
            (calc.set_positions must have been called before this point!).
        wfs: WaveFunctions
            Class taking care of wave-functions, projectors, k-point related
            quantities and symmetries.
        poisson_solver: PoissonSolver
            Multigrid or FFT poisson solver (not required if the
            ``Perturbation`` to be solved for has a ``solve_poisson`` member
            function). 
        dtype: ...
            dtype of the density response.
            
        """
        
        # Store ground-state quantities
        self.hamiltonian = calc.hamiltonian
        self.density = calc.density

        self.wfs = wfs

        # Get list of k-point containers
        self.kpt_u = wfs.kpt_u

        # Poisson solver
        if poisson_solver is None:
            # Solver must be provided by the perturbation
            self.poisson = None
            self.solve_poisson = None
        else:
            self.poisson = poisson_solver
            self.solve_poisson = self.poisson.solve_neutral
       
        # Store grid-descriptors
        self.gd = calc.density.gd
        self.finegd = calc.density.finegd

        # dtype for ground-state wave-functions
        self.gs_dtype = calc.wfs.dtype
        # dtype for the perturbing potential and density
        self.dtype = dtype
        
        # Grid transformer -- convert array from coarse to fine grid
        self.interpolator = Transformer(self.gd, self.finegd, nn=3,
                                        dtype=self.dtype)
        # Grid transformer -- convert array from fine to coarse grid
        self.restrictor = Transformer(self.finegd, self.gd, nn=3,
                                      dtype=self.dtype)

        # Sternheimer operator
        self.sternheimer_operator = None
        # Krylov solver
        self.linear_solver = None

        # Phases for transformer objects - since the perturbation determines
        # the form of the density response this is obtained from the
        # perturbation in the ``__call__`` member function below.
        self.phase_cd = None

        # Array attributes
        self.nt1_G = None
        self.vHXC1_G = None        
        self.nt1_g = None
        self.vH1_g = None

        # Perturbation
        self.perturbation = None
        
        # Number of occupied bands
        nvalence = calc.wfs.nvalence
        self.nbands = nvalence / 2 + nvalence % 2
        assert self.nbands <= calc.wfs.bd.nbands
                                  
        self.initialized = False

        self.parameters = {}
        self.set(**kwargs)

    def clean(self):
        """Cleanup before call to ``__call__``."""

        self.perturbation = None
        self.solve_poisson = None

        self.nt1_G = None
        self.vHXC1_G = None        
        self.nt1_g = None
        self.vH1_g = None
        
    def __call__(self, perturbation):
        """Calculate density response (derivative) to perturbation.

        Parameters
        ----------
        perturbation: Perturbation
            Class implementing the perturbing potential. Must provide an
            ``apply`` member function implementing the multiplication of the
            perturbing potential to a (set of) state vector(s).
            
        """
        
        assert self.initialized, ("Linear response calculator "
                                  "not initizalized.")
        self.clean()
        
        if self.poisson is None:
            assert hasattr(perturbation, 'solve_poisson')
            self.solve_poisson = perturbation.solve_poisson

        # Store perturbation - used in other member functions
        self.perturbation = perturbation
        # Reset mixer
        self.mixer.reset()
        # Reset wave-functions
        self.wfs.reset()

        # Set phase attribute for Transformer objects
        self.phase_cd = self.perturbation.get_phase_cd()

        # Parameters for the SC-loop
        p = self.parameters
        max_iter = p['max_iter']
        tolerance = p['tolerance_sc']
        
        for iter in range(max_iter):

            if iter == 0:
                self.first_iteration()
            else:
                print("iter:%3i\t" % iter, end=' ')
                norm = self.iteration()
                print("abs-norm: %6.3e\t" % norm, end=' ')
                print(("integrated density response (abs): % 5.2e (%5.2e) "
                       % (self.gd.integrate(self.nt1_G.real),
                          self.gd.integrate(np.absolute(self.nt1_G)))))
                       
                if norm < tolerance:
                    print(("self-consistent loop converged in %i iterations"
                           % iter))
                    break
                
            if iter == max_iter:
                raise RuntimeError, ("self-consistent loop did not converge "
                                     "in %i iterations" % iter)
   
    def set(self, **kwargs):
        """Set parameters for calculation."""

        # Check for legal input parameters
        for key, value in kwargs.items():
            if not key in ResponseCalculator.parameters:
                raise TypeError("Unknown keyword argument: '%s'" % key)

        # Insert default values if not given
        for key, value in ResponseCalculator.parameters.items():
            if key not in kwargs:
                kwargs[key] = value

        self.parameters.update(kwargs)
            
    def initialize(self, spos_ac):
        """Make the object ready for a calculation."""

        # Parameters
        p = self.parameters
        beta = p['beta']
        nmaxold = p['nmaxold']
        weight = p['weight']
        use_pc = p['use_pc']
        tolerance_sternheimer = p['tolerance_sternheimer']
        max_iter_krylov = p['max_iter_krylov']
        krylov_solver = p['krylov_solver']
                
        # Initialize WaveFunctions attribute
        self.wfs.initialize(spos_ac)
        
        # Initialize mixer
        # weight = 1 -> no metric is used
        self.mixer = BaseMixer(beta=beta, nmaxold=nmaxold,
                               weight=weight, dtype=self.dtype)
        self.mixer.initialize_metric(self.gd)
        
        # Linear operator in the Sternheimer equation
        self.sternheimer_operator = \
            SternheimerOperator(self.hamiltonian, self.wfs, self.gd,
                                dtype=self.gs_dtype)

        # Preconditioner for the Sternheimer equation
        if p['use_pc']:
            pc = ScipyPreconditioner(self.gd,
                                     self.sternheimer_operator.project,
                                     dtype=self.gs_dtype)
        else:
            pc = None

        #XXX K-point of the pc must be set in the k-point loop -> store a ref.
        self.pc = pc
        # Linear solver for the solution of Sternheimer equation            
        self.linear_solver = ScipyLinearSolver(method=krylov_solver,
                                               preconditioner=pc,
                                               tolerance=tolerance_sternheimer,
                                               max_iter=max_iter_krylov)

        self.initialized = True

    def first_iteration(self):
        """Perform first iteration of sc-loop."""

        self.wave_function_variations()
        self.density_response()
        self.mixer.mix(self.nt1_G, [], phase_cd=self.phase_cd)
        self.interpolate_density()
        
    def iteration(self):
        """Perform iteration."""

        # Update variation in the effective potential
        self.effective_potential_variation()
        # Update wave function variations
        self.wave_function_variations()
        # Update density
        self.density_response()
        # Mix - supply phase_cd here for metric inside the mixer
        self.mixer.mix(self.nt1_G, [], phase_cd=self.phase_cd)
        norm = self.mixer.get_charge_sloshing()

        self.interpolate_density()
       
        return norm

    def interpolate_density(self):
        """Interpolate density derivative onto the fine grid."""

        self.nt1_g = self.finegd.zeros(dtype=self.dtype)
        self.interpolator.apply(self.nt1_G, self.nt1_g, phases=self.phase_cd)
        
    def effective_potential_variation(self):
        """Calculate derivative of the effective potential (Hartree + XC)."""

        # Hartree part
        vHXC1_g = self.finegd.zeros(dtype=self.dtype)
        self.solve_poisson(vHXC1_g, self.nt1_g)
        # Store for evaluation of second order derivative
        self.vH1_g = vHXC1_g.copy()
        
        # XC part
        nt_sg = self.density.nt_sg
        fxct_sg = np.zeros_like(nt_sg)
        self.hamiltonian.xc.calculate_fxc(self.finegd, nt_sg, fxct_sg)
        vHXC1_g += fxct_sg[0] * self.nt1_g

        # Transfer to coarse grid
        self.vHXC1_G = self.gd.zeros(dtype=self.dtype)
        self.restrictor.apply(vHXC1_g, self.vHXC1_G, phases=self.phase_cd)
    
    def wave_function_variations(self):
        """Calculate variation in the wave-functions.

        Parameters
        ----------
        v1_G: ndarray
            Variation of the local effective potential (Hartree + XC).

        """

        verbose = self.parameters['verbose']

        if verbose:
            print("Calculating wave function variations")

        if self.perturbation.has_q():
            q_c = self.perturbation.get_q()
            kplusq_k = self.wfs.kd.find_k_plus_q(q_c)
        else:
            kplusq_k = None

        # Calculate wave-function variations for all k-points.
        for kpt in self.kpt_u:

            k = kpt.k

            if verbose:
                print("k-point %2.1i" % k)
            
            # Index of k+q vector
            if kplusq_k is None:
                kplusq = k
                kplusqpt = kpt
            else:
                kplusq = kplusq_k[k]
                kplusqpt = self.kpt_u[kplusq]

            # Ground-state and first-order wave-functions
            psit_nG = kpt.psit_nG
            psit1_nG = kpt.psit1_nG
            # Update the SternheimerOperator
            self.sternheimer_operator.set_k(k)
            self.sternheimer_operator.set_kplusq(kplusq)
            # Update preconditioner
            if self.pc is not None:
                # k+q
                self.pc.set_kpt(kplusqpt)
                
            # Right-hand side of Sternheimer equations
            # k and k+q
            # XXX should only be done once for all k-points but maybe too cheap
            # to bother ??
            rhs_nG = self.gd.zeros(n=self.nbands, dtype=self.gs_dtype)            
            self.perturbation.apply(psit_nG, rhs_nG, self.wfs, k, kplusq)
            if self.vHXC1_G is not None:
                rhs_nG += self.vHXC1_G * psit_nG
            # Project out occupied subspace
            self.sternheimer_operator.project(rhs_nG)
              
            # Loop over occupied bands
            for n in range(self.nbands):

                # Update band index in SternheimerOperator
                self.sternheimer_operator.set_band(n)
                # Get view of the Bloch function derivative
                psit1_G = psit1_nG[n]
                # Rhs of Sternheimer equation                
                rhs_G = -1 * rhs_nG[n]
               
                # Solve Sternheimer equation
                iter, info = self.linear_solver.solve(self.sternheimer_operator,
                                                      psit1_G, rhs_G)
                
                if verbose:
                    print("\tBand %2.1i -" % n, end=' ')
                    
                if info == 0:
                    if verbose:
                        print("linear solver converged in %i iterations" % iter)
                elif info > 0:
                    assert False, ("linear solver did not converge in maximum "
                                   "number (=%i) of iterations for "
                                   "k-point number %d" % (iter, k))
                else:
                    assert False, ("linear solver failed to converge")

    def density_response(self):
        """Calculate density response from variation in the wave-functions."""

        # Density might be complex
        self.nt1_G = self.gd.zeros(dtype=self.dtype)

        for kpt in self.kpt_u:
            # The weight of the k-points includes spin-degeneracy
            w = kpt.weight
            # Wave functions
            psit_nG = kpt.psit_nG
            psit1_nG = kpt.psit1_nG

            for psit_G, psit1_G in zip(psit_nG, psit1_nG):
                # NOTICE: this relies on the automatic down-cast of the complex
                # array on the rhs to a real array when the lhs is real !!
                # Factor 2 for time-reversal symmetry
                self.nt1_G += 2 * w * psit_G.conj() * psit1_G
