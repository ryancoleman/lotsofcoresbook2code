# Written by Lauri Lehtovaara, 2007

"""This module implements time propagators for time-dependent density
functional theory calculations."""

import sys

import numpy as np

from gpaw.utilities.blas import axpy
from gpaw.utilities.blas import dotc

from gpaw.mpi import rank, run

from gpaw.tddft.utils import MultiBlas
from gpaw.tddft.tdopers import DummyDensity

###############################################################################
# DummyKPoint
###############################################################################
class DummyKPoint(object):
    __slots__ = ('psit_nG',)

###############################################################################
# DummyPropagator
###############################################################################
class DummyPropagator:
    """Time propagator
    
    The DummyPropagator-class is the VIRTUAL base class for all propagators.
    
    """
    def __init__(self, td_density, td_hamiltonian, td_overlap,
                solver, preconditioner, gd, timer):
        """Create the DummyPropagator-object.
        
        Parameters
        ----------
        td_density: TimeDependentDensity
            the time-dependent density
        td_hamiltonian: TimeDependentHamiltonian
            the time-dependent hamiltonian
        td_overlap: TimeDependentOverlap
            the time-dependent overlap operator
        solver: LinearSolver
            solver for linear equations
        preconditioner: Preconditioner
            preconditioner for linear equations
        gd: GridDescriptor
            coarse (/wavefunction) grid descriptor
        timer: Timer
            timer
        
        """
        self.td_density = td_density
        self.td_hamiltonian = td_hamiltonian
        self.td_overlap = td_overlap

        self.wfs = td_density.get_wavefunctions()

        self.solver = solver
        self.preconditioner = preconditioner
        self.gd = gd
        self.timer = timer

        self.mblas = MultiBlas(gd)

    # Solve M psin = psi
    def apply_preconditioner(self, psi, psin):
        """Solves preconditioner equation.
        
        Parameters
        ----------
        psi: List of coarse grids
            the known wavefunctions
        psin: List of coarse grids
            the result
        
        """
        self.timer.start('Solve TDDFT preconditioner')
        if self.preconditioner is not None:
            self.preconditioner.apply(self.kpt, psi, psin)
        else:
            psin[:] = psi
        self.timer.stop('Solve TDDFT preconditioner')


    def propagate(self, time, time_step):
        """Propagate wavefunctions once. 
        
        Parameters
        ----------
        time: float
            the current time
        time_step: float
            the time step
        
        """
        raise RuntimeError('Error in DummyPropagator: '
                            'Member function propagate is virtual.')


###############################################################################
# ExplicitCrankNicolson
###############################################################################
class ExplicitCrankNicolson(DummyPropagator):
    """Explicit Crank-Nicolson propagator
    
    Crank-Nicolson propagator, which approximates the time-dependent
    Hamiltonian to be unchanged during one iteration step.
    
    (S(t) + .5j dt H(t) / hbar) psi(t+dt) = (S(t) - .5j dt H(t) / hbar) psi(t)
    
    """
    
    def __init__(self, td_density, td_hamiltonian, td_overlap,
                 solver, preconditioner, gd, timer):
        """Create ExplicitCrankNicolson-object.
        
        Parameters
        ----------
        td_density: TimeDependentDensity
            time-dependent density
        td_hamiltonian: TimeDependentHamiltonian
            time-dependent hamiltonian
        td_overlap: TimeDependentOverlap
            time-dependent overlap operator
        solver: LinearSolver
            solver for linear equations
        preconditioner: Preconditioner
            preconditioner for linear equations
        gd: GridDescriptor
            coarse (/wavefunction) grid descriptor
        timer: Timer
            timer
        
        """
        DummyPropagator.__init__(self, td_density, td_hamiltonian, td_overlap,
                                 solver, preconditioner, gd, timer)

        self.tmp_kpt_u = None
        self.hpsit = None
        self.spsit = None
        self.sinvhpsit = None
        
    # ( S + i H dt/2 ) psit(t+dt) = ( S - i H dt/2 ) psit(t)
    def propagate(self, time, time_step):
        """Propagate wavefunctions. 
        
        Parameters
        ----------
        time: float
            the current time
        time_step: float
            time step
        
        """
        self.niter = 0

        # Allocate temporary wavefunctions
        if self.tmp_kpt_u is None:
            self.tmp_kpt_u = []
            for kpt in self.wfs.kpt_u:
                tmp_kpt = DummyKPoint()
                tmp_kpt.psit_nG = self.gd.empty(n=len(kpt.psit_nG),
                                                dtype=complex)
                self.tmp_kpt_u.append(tmp_kpt)

        # Allocate memory for Crank-Nicolson stuff
        nvec = len(self.wfs.kpt_u[0].psit_nG)
        if self.hpsit is None:
            self.hpsit = self.gd.zeros(nvec, dtype=complex)
        if self.spsit is None:
            self.spsit = self.gd.zeros(nvec, dtype=complex)


        self.timer.start('Update time-dependent operators')

        # Update overlap S(t) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t) based on the wavefunctions psit_nG
        # in kpt_u for t = time. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Update Hamiltonian H(t) to reflect density rho(t)
        self.td_hamiltonian.update(self.td_density.get_density(), time)

        self.timer.stop('Update time-dependent operators')

        # Copy current wavefunctions psit_nG to work wavefunction arrays
        for u, kpt in enumerate(self.wfs.kpt_u):
            self.tmp_kpt_u[u].psit_nG[:] = kpt.psit_nG


        # Euler step
        # Overwrite psit_nG in tmp_kpt_u by (1 - i S^(-1)(t) H(t) dt) psit_nG
        # from corresponding kpt_u in a Euler step before predicting psit(t+dt)
        for [kpt, rhs_kpt] in zip(self.wfs.kpt_u, self.tmp_kpt_u):
            self.solve_propagation_equation(kpt, rhs_kpt, time_step, guess=True)

        # update projections before exiting
        self.td_overlap.update(self.wfs)

        return self.niter


    # ( S + i H dt/2 ) psit(t+dt) = ( S - i H dt/2 ) psit(t)
    def solve_propagation_equation(self, kpt, rhs_kpt, time_step, guess=False):

        # kpt is guess, rhs_kpt is used to calculate rhs and is overwritten
        nvec = len(rhs_kpt.psit_nG)

        assert kpt != rhs_kpt, 'Data race condition detected'
        assert len(kpt.psit_nG) == nvec, 'Incompatible lhs/rhs vectors'

        self.timer.start('Apply time-dependent operators')
        # Store H psi(t) as hpsit and S psit(t) as spsit
        self.td_overlap.update_k_point_projections(self.wfs, kpt, rhs_kpt.psit_nG)
        self.td_hamiltonian.apply(kpt, rhs_kpt.psit_nG, self.hpsit,
                                  calculate_P_ani=False)
        self.td_overlap.apply(rhs_kpt.psit_nG, self.spsit, self.wfs, kpt,
                              calculate_P_ani=False)
        self.timer.stop('Apply time-dependent operators')

        # Update rhs_kpt.psit_nG to reflect ( S - i H dt/2 ) psit(t)
        #rhs_kpt.psit_nG[:] = self.spsit - .5J * self.hpsit * time_step
        rhs_kpt.psit_nG[:] = self.spsit
        self.mblas.multi_zaxpy(-.5j*time_step, self.hpsit, rhs_kpt.psit_nG, nvec)

        if guess:
            if self.sinvhpsit is None:
                self.sinvhpsit = self.gd.zeros(len(kpt.psit_nG), dtype=complex)

            # Update estimate of psit(t+dt) to ( 1 - i S^(-1) H dt ) psit(t)
            self.td_overlap.apply_inverse(self.hpsit, self.sinvhpsit, self.wfs, kpt, use_cg=False)
            self.mblas.multi_zaxpy(-1.0j*time_step, self.sinvhpsit,
                                   kpt.psit_nG, nvec)

        # Information needed by solver.solve -> self.dot
        self.kpt = kpt
        self.time_step = time_step

        # Solve A x = b where A is (S + i H dt/2) and b = rhs_kpt.psit_nG
        self.niter += self.solver.solve(self, kpt.psit_nG, rhs_kpt.psit_nG)

    # ( S + i H dt/2 ) psi
    def dot(self, psi, psin):
        """Applies the propagator matrix to the given wavefunctions.

        Parameters
        ----------
        psi: List of coarse grids
            the known wavefunctions
        psin: List of coarse grids
            the result ( S + i H dt/2 ) psi

        """
        self.timer.start('Apply time-dependent operators')
        self.td_overlap.update_k_point_projections(self.wfs, self.kpt, psi)
        self.td_hamiltonian.apply(self.kpt, psi, self.hpsit,
                                  calculate_P_ani=False)
        self.td_overlap.apply(psi, self.spsit, self.wfs, self.kpt, calculate_P_ani=False)
        self.timer.stop('Apply time-dependent operators')

        # psin[:] = self.spsit + .5J * self.time_step * self.hpsit
        psin[:] = self.spsit
        self.mblas.multi_zaxpy(.5j*self.time_step, self.hpsit, psin, len(psi))



###############################################################################
# SemiImplicitCrankNicolson
###############################################################################
class SemiImplicitCrankNicolson(ExplicitCrankNicolson):
    """Semi-implicit Crank-Nicolson propagator
    
    Crank-Nicolson propagator, which first approximates the time-dependent
    Hamiltonian to be unchanged during one iteration step to predict future
    wavefunctions. Then the approximations for the future wavefunctions are
    used to approximate the Hamiltonian at the middle of the time step.
    
    (S(t) + .5j dt H(t) / hbar) psi(t+dt) = (S(t) - .5j dt H(t) / hbar) psi(t)
    (S(t) + .5j dt H(t+dt/2) / hbar) psi(t+dt) 
    = (S(t) - .5j dt H(t+dt/2) / hbar) psi(t)
    
    """
    
    def __init__(self, td_density, td_hamiltonian, td_overlap, 
                 solver, preconditioner, gd, timer):
        """Create SemiImplicitCrankNicolson-object.
        
        Parameters
        ----------
        td_density: TimeDependentDensity
            the time-dependent density
        td_hamiltonian: TimeDependentHamiltonian
            the time-dependent hamiltonian
        td_overlap: TimeDependentOverlap
            the time-dependent overlap operator
        solver: LinearSolver
            solver for linear equations
        preconditioner: Preconditioner
            preconditioner for linear equations
        gd: GridDescriptor
            coarse (wavefunction) grid descriptor
        timer: Timer
            timer
        
        """
        ExplicitCrankNicolson.__init__(self, td_density, td_hamiltonian,
                          td_overlap, solver, preconditioner, gd, timer)

        self.old_kpt_u = None


    def propagate(self, time, time_step):
        """Propagate wavefunctions once.
        
        Parameters
        ----------
        time: float
            the current time
        time_step: float
            time step
        """

        self.niter = 0

        # Allocate old/temporary wavefunctions
        if self.old_kpt_u is None:
            self.old_kpt_u = []
            for kpt in self.wfs.kpt_u:
                old_kpt = DummyKPoint()
                old_kpt.psit_nG = self.gd.empty(n=len(kpt.psit_nG),
                                                dtype=complex)
                self.old_kpt_u.append(old_kpt)

        if self.tmp_kpt_u is None:
            self.tmp_kpt_u = []
            for kpt in self.wfs.kpt_u:
                tmp_kpt = DummyKPoint()
                tmp_kpt.psit_nG = self.gd.empty(n=len(kpt.psit_nG),
                                                dtype=complex)
                self.tmp_kpt_u.append(tmp_kpt)

        # Allocate memory for Crank-Nicolson stuff
        nvec = len(self.wfs.kpt_u[0].psit_nG)
        if self.hpsit is None:
            self.hpsit = self.gd.zeros(nvec, dtype=complex)
        if self.spsit is None:
            self.spsit = self.gd.zeros(nvec, dtype=complex)


        self.timer.start('Update time-dependent operators')

        # Update overlap S(t) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t) based on the wavefunctions psit_nG
        # in kpt_u for t = time. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Update Hamiltonian H(t) to reflect density rho(t)
        self.td_hamiltonian.update(self.td_density.get_density(), time)

        self.timer.stop('Update time-dependent operators')

        # Copy current wavefunctions psit_nG to work and old wavefunction arrays
        for u, kpt in enumerate(self.wfs.kpt_u):
            self.old_kpt_u[u].psit_nG[:] = kpt.psit_nG
            self.tmp_kpt_u[u].psit_nG[:] = kpt.psit_nG


        # Predictor step
        # Overwrite psit_nG in tmp_kpt_u by (1 - i S^(-1)(t) H(t) dt) psit_nG
        # from corresponding kpt_u in a Euler step before predicting psit(t+dt)
        for [kpt, rhs_kpt] in zip(self.wfs.kpt_u, self.tmp_kpt_u):
            self.solve_propagation_equation(kpt, rhs_kpt, time_step, guess=True)

        self.timer.start('Update time-dependent operators')

        # Update overlap S(t+dt) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t+dt) based on the wavefunctions psit_nG in
        # kpt_u for t = time+time_step. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Estimate Hamiltonian H(t+dt/2) by averaging H(t) and H(t+dt)
        # and retain the difference for a half-way Hamiltonian dH(t+dt/2).
        self.td_hamiltonian.half_update(self.td_density.get_density(),
                                        time + time_step)

        # Estimate overlap S(t+dt/2) by averaging S(t) and S(t+dt) #TODO!!!
        self.td_overlap.half_update(self.wfs)

        self.timer.stop('Update time-dependent operators')

        # Corrector step
        # Use predicted psit_nG in kpt_u as an initial guess, whereas the old 
        # wavefunction in old_kpt_u are used to calculate rhs based on psit(t)
        for [kpt, rhs_kpt] in zip(self.wfs.kpt_u, self.old_kpt_u):
            # Average of psit(t) and predicted psit(t+dt)
            mean_psit_nG = 0.5*(kpt.psit_nG + rhs_kpt.psit_nG)
            self.td_hamiltonian.half_apply(kpt, mean_psit_nG, self.hpsit)
            self.td_overlap.apply_inverse(self.hpsit, self.sinvhpsit, self.wfs, kpt, use_cg=False)

            # Update kpt.psit_nG to reflect psit(t+dt) - i S^(-1) dH(t+dt/2) dt/2 psit(t+dt/2)
            kpt.psit_nG[:] = kpt.psit_nG - .5J * self.sinvhpsit * time_step
            self.mblas.multi_zaxpy(-.5j*time_step, self.sinvhpsit, kpt.psit_nG, nvec)

            self.solve_propagation_equation(kpt, rhs_kpt, time_step)
 
        # update projections before exiting
        self.td_overlap.update(self.wfs)
       
        return self.niter


    # ( S + i H dt/2 ) psit(t+dt) = ( S - i H dt/2 ) psit(t)
    def solve_propagation_equation(self, kpt, rhs_kpt, time_step, guess=False):

        # kpt is guess, rhs_kpt is used to calculate rhs and is overwritten
        nvec = len(rhs_kpt.psit_nG)

        assert kpt != rhs_kpt, 'Data race condition detected'
        assert len(kpt.psit_nG) == nvec, 'Incompatible lhs/rhs vectors'

        self.timer.start('Apply time-dependent operators')
        # Store H psi(t) as hpsit and S psit(t) as spsit
        self.td_overlap.update_k_point_projections(self.wfs, kpt, rhs_kpt.psit_nG)
        self.td_hamiltonian.apply(kpt, rhs_kpt.psit_nG, self.hpsit,
                                  calculate_P_ani=False)
        self.td_overlap.apply(rhs_kpt.psit_nG, self.spsit, self.wfs, kpt,
                              calculate_P_ani=False)
        self.timer.stop('Apply time-dependent operators')

        #self.mblas.multi_zdotc(self.shift, rhs_kpt.psit_nG, self.hpsit, nvec)
        #self.shift *= self.gd.dv
        #self.mblas.multi_zdotc(self.tmp_shift, rhs_kpt.psit_nG, self.spsit, nvec)
        #self.tmp_shift *= self.gd.dv
        #self.shift /= self.tmp_shift

        # Update rhs_kpt.psit_nG to reflect ( S - i H dt/2 ) psit(t)
        #rhs_kpt.psit_nG[:] = self.spsit - .5J * self.hpsit * time_step
        rhs_kpt.psit_nG[:] = self.spsit
        self.mblas.multi_zaxpy(-.5j*time_step, self.hpsit, rhs_kpt.psit_nG, nvec)
        # Apply shift -i eps S t/2
        #self.mblas.multi_zaxpy(-.5j*time_step * (-self.shift), self.spsit, rhs_kpt.psit_nG, nvec)

        if guess:
            if self.sinvhpsit is None:
                self.sinvhpsit = self.gd.zeros(len(kpt.psit_nG), dtype=complex)

            # Update estimate of psit(t+dt) to ( 1 - i S^(-1) H dt ) psit(t)
            self.td_overlap.apply_inverse(self.hpsit, self.sinvhpsit, self.wfs, kpt, use_cg=False)
            self.mblas.multi_zaxpy(-1.0j*time_step, self.sinvhpsit,
                                   kpt.psit_nG, nvec)

        # Information needed by solver.solve -> self.dot
        self.kpt = kpt
        self.time_step = time_step

        # Solve A x = b where A is (S + i H dt/2) and b = rhs_kpt.psit_nG
        self.niter += self.solver.solve(self, kpt.psit_nG, rhs_kpt.psit_nG)

        # Apply shift exp(i eps t)
        #self.phase_shift = np.exp(1.0J * self.shift * time_step)
        #self.mblas.multi_scale(self.phase_shift, kpt.psit_nG, nvec)

    # ( S + i H dt/2 ) psi
    def dot(self, psi, psin):
        """Applies the propagator matrix to the given wavefunctions.

        Parameters
        ----------
        psi: List of coarse grids
            the known wavefunctions
        psin: List of coarse grids
            the result ( S + i H dt/2 ) psi

        """
        ExplicitCrankNicolson.dot(self, psi, psin)
        # Apply shift -i eps S t/2 
        #self.mblas.multi_zaxpy(.5j * self.time_step * (-self.shift), self.spsit, psin, len(psi))

class EhrenfestPAWSICN(ExplicitCrankNicolson):
    """Semi-implicit Crank-Nicolson propagator for Ehrenfest dynamics
       TODO: merge this with the ordinary SICN
    """
    def __init__(self, td_density, td_hamiltonian, td_overlap, 
                 solver, preconditioner, gd, timer, corrector_guess = True,
                 predictor_guess = (True,False), use_cg = (False,False)):
        """Create SemiImplicitCrankNicolson-object.
        
        Parameters
        ----------
        td_density: TimeDependentDensity
            the time-dependent density
        td_hamiltonian: TimeDependentHamiltonian
            the time-dependent hamiltonian
        td_overlap: TimeDependentOverlap
            the time-dependent overlap operator
        solver: LinearSolver
            solver for linear equations
        preconditioner: Preconditioner
            preconditioner for linear equations
        gd: GridDescriptor
            coarse (wavefunction) grid descriptor
        timer: Timer
            timer
        corrector_guess: Bool
            use initial guess for the corrector step (default is True)
        predictor_guess: (Bool, Bool)
            use (first, second) order initial guesses for the predictor step
            default is (True, False)
        use_cg: (Bool, Bool)
            use CG for calculating the inverse overlap (predictor, corrector)
            default is (False, False)
        
        """
        ExplicitCrankNicolson.__init__(self, td_density, td_hamiltonian,
                          td_overlap, solver, preconditioner, gd, timer)

        self.old_kpt_u = None
        self.corrector_guess = corrector_guess
        self.predictor_guess = predictor_guess
        self.use_cg = use_cg

        #self.hsinvhpsit = None
        self.sinvh2psit = None

    def update_velocities(self, v_at_new, v_at_old = None):
        self.v_at = v_at_new.copy()
        if(v_at_old is not None):
            self.v_at_old = v_at_old.copy()

    def propagate(self, time, time_step, v_a):
        """Propagate wavefunctions once.
        
        Parameters
        ----------
        time: float
            the current time
        time_step: float
            time step
        """

        self.niter = 0

        #update the atomic velocities which are required
        #for calculating the P term
        self.update_velocities(v_a)

        # Allocate old/temporary wavefunctions
        if self.old_kpt_u is None:
            self.old_kpt_u = []
            for kpt in self.wfs.kpt_u:
                old_kpt = DummyKPoint()
                old_kpt.psit_nG = self.gd.empty(n=len(kpt.psit_nG),
                                                dtype=complex)
                self.old_kpt_u.append(old_kpt)

        if self.tmp_kpt_u is None:
            self.tmp_kpt_u = []
            for kpt in self.wfs.kpt_u:
                tmp_kpt = DummyKPoint()
                tmp_kpt.psit_nG = self.gd.empty(n=len(kpt.psit_nG),
                                                dtype=complex)
                self.tmp_kpt_u.append(tmp_kpt)

        # Allocate memory for Crank-Nicolson stuff
        nvec = len(self.wfs.kpt_u[0].psit_nG)
        if self.hpsit is None:
            self.hpsit = self.gd.zeros(nvec, dtype=complex)
        if self.spsit is None:
            self.spsit = self.gd.zeros(nvec, dtype=complex)


        self.timer.start('Update time-dependent operators')

        # Update overlap S(t) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t) based on the wavefunctions psit_nG
        # in kpt_u for t = time. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Update Hamiltonian H(t) to reflect density rho(t)
        self.td_hamiltonian.update(self.td_density.get_density(), time)

        self.timer.stop('Update time-dependent operators')

        # Copy current wavefunctions psit_nG to work and old wavefunction arrays
        for u, kpt in enumerate(self.wfs.kpt_u):
            self.old_kpt_u[u].psit_nG[:] = kpt.psit_nG
            self.tmp_kpt_u[u].psit_nG[:] = kpt.psit_nG

        #print 'P_ani[0] =', self.wfs.kpt_u[0].P_ani[0]
        #print self.test


        # Predictor step
        # Overwrite psit_nG in tmp_kpt_u by (1 - i S^(-1)(t) H(t) dt) psit_nG
        # from corresponding kpt_u in a Euler step before predicting psit(t+dt)
        #self.v_at = self.v_at_old.copy() #v(t) for predictor step
        for [kpt, rhs_kpt] in zip(self.wfs.kpt_u, self.tmp_kpt_u):
            #print 'self.predictor_guess[0]', self.predictor_guess[0]
            self.solve_propagation_equation(kpt, rhs_kpt, time_step, guess=self.predictor_guess[0])

        self.timer.start('Update time-dependent operators')

        # Update overlap S(t+dt) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t+dt) based on the wavefunctions psit_nG in
        # kpt_u for t = time+time_step. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Estimate Hamiltonian H(t+dt/2) by averaging H(t) and H(t+dt)
        # and retain the difference for a half-way Hamiltonian dH(t+dt/2).
        self.td_hamiltonian.half_update(self.td_density.get_density(),
                                        time + time_step)

        # Estimate overlap S(t+dt/2) by averaging S(t) and S(t+dt) #TODO!!!
        self.td_overlap.half_update(self.wfs)

        self.timer.stop('Update time-dependent operators')

        # Corrector step
        # Use predicted psit_nG in kpt_u as an initial guess, whereas the old 
        # wavefunction in old_kpt_u are used to calculate rhs based on psit(t)
        for [kpt, rhs_kpt] in zip(self.wfs.kpt_u, self.old_kpt_u):
            # Average of psit(t) and predicted psit(t+dt)
            if(self.corrector_guess):
                mean_psit_nG = 0.5*(kpt.psit_nG + rhs_kpt.psit_nG)
                self.td_hamiltonian.half_apply(kpt, mean_psit_nG, self.hpsit)
                self.td_overlap.apply_inverse(self.hpsit, self.sinvhpsit, self.wfs, kpt, use_cg=self.use_cg[1])

                # Update kpt.psit_nG to reflect psit(t+dt) - i S^(-1) dH(t+dt/2) dt/2 psit(t+dt/2)
                kpt.psit_nG[:] = kpt.psit_nG - .5J * self.sinvhpsit * time_step
                self.mblas.multi_zaxpy(-.5j*time_step, self.sinvhpsit, kpt.psit_nG, nvec)

            self.solve_propagation_equation(kpt, rhs_kpt, time_step)
 
        # update projections before exiting
        self.td_overlap.update(self.wfs)
       
        return self.niter


    # ( S + i H dt/2 ) psit(t+dt) = ( S - i H dt/2 ) psit(t)
    def solve_propagation_equation(self, kpt, rhs_kpt, time_step, calculate_P_ani=False, guess=False):

        # kpt is guess, rhs_kpt is used to calculate rhs and is overwritten
        nvec = len(rhs_kpt.psit_nG)

        assert kpt != rhs_kpt, 'Data race condition detected'
        assert len(kpt.psit_nG) == nvec, 'Incompatible lhs/rhs vectors'

        self.timer.start('Apply time-dependent operators')
        # Store H psi(t) as hpsit and S psit(t) as spsit
        self.td_overlap.update_k_point_projections(self.wfs, kpt, rhs_kpt.psit_nG)
        self.td_hamiltonian.apply(kpt, rhs_kpt.psit_nG, self.hpsit,
                                  calculate_P_ani=False)
        self.td_hamiltonian.calculate_paw_correction(rhs_kpt.psit_nG, self.hpsit, self.wfs, kpt, self.v_at, calculate_P_ani=False)
            
        self.td_overlap.apply(rhs_kpt.psit_nG, self.spsit, self.wfs, kpt,
                              calculate_P_ani=False)
        self.timer.stop('Apply time-dependent operators')

        #self.mblas.multi_zdotc(self.shift, rhs_kpt.psit_nG, self.hpsit, nvec)
        #self.shift *= self.gd.dv
        #self.mblas.multi_zdotc(self.tmp_shift, rhs_kpt.psit_nG, self.spsit, nvec)
        #self.tmp_shift *= self.gd.dv
        #self.shift /= self.tmp_shift

        # Update rhs_kpt.psit_nG to reflect ( S - i H dt/2 ) psit(t)
        #rhs_kpt.psit_nG[:] = self.spsit - .5J * self.hpsit * time_step
        rhs_kpt.psit_nG[:] = self.spsit
        self.mblas.multi_zaxpy(-.5j*time_step, self.hpsit, rhs_kpt.psit_nG, nvec)
        # Apply shift -i eps S t/2
        #self.mblas.multi_zaxpy(-.5j*time_step * (-self.shift), self.spsit, rhs_kpt.psit_nG, nvec)

        if guess:
            if self.sinvhpsit is None:
                self.sinvhpsit = self.gd.zeros(len(kpt.psit_nG), dtype=complex)

            if self.predictor_guess[1]:
                if self.sinvh2psit is None:
                    self.sinvh2psit = self.gd.zeros(len(kpt.psit_nG), dtype=complex)

            # Update estimate of psit(t+dt) to ( 1 - i S^(-1) H dt ) psit(t)
            #print 'self.use_cg[0]', self.use_cg[0]
            self.td_overlap.apply_inverse(self.hpsit, self.sinvhpsit, self.wfs, kpt, use_cg=self.use_cg[0])
            #assert self.use_cg[0] is False
            #self.td_overlap.apply_inverse(self.hpsit, self.sinvhpsit, self.wfs, kpt, use_cg=False)
            
            self.mblas.multi_zaxpy(-1.0j*time_step, self.sinvhpsit,
                                   kpt.psit_nG, nvec)
            #print 'using guess for P step'
            if(self.predictor_guess[1]):
                #print 'using 2nd order guess for P step'
                self.td_hamiltonian.apply(kpt, self.sinvhpsit, self.sinvh2psit,
                                  calculate_P_ani=False)
                self.td_hamiltonian.calculate_paw_correction(self.sinvhpsit, self.sinvh2psit, self.wfs, kpt, self.v_at, calculate_P_ani=False)
                self.td_overlap.apply_inverse(self.sinvh2psit, self.sinvh2psit, self.wfs, kpt, use_cg=self.use_cg[0])
                self.mblas.multi_zaxpy(-.5 * time_step*time_step, self.sinvh2psit,
                                   kpt.psit_nG, nvec)
                

        # Information needed by solver.solve -> self.dot
        self.kpt = kpt
        self.time_step = time_step

        # Solve A x = b where A is (S + i H dt/2) and b = rhs_kpt.psit_nG
        self.niter += self.solver.solve(self, kpt.psit_nG, rhs_kpt.psit_nG)

        # Apply shift exp(i eps t)
        #self.phase_shift = np.exp(1.0J * self.shift * time_step)
        #self.mblas.multi_scale(self.phase_shift, kpt.psit_nG, nvec)

    # ( S + i H dt/2 ) psi
    def dot(self, psi, psin):
        """Applies the propagator matrix to the given wavefunctions.

        Parameters
        ----------
        psi: List of coarse grids
            the known wavefunctions
        psin: List of coarse grids
            the result ( S + i H dt/2 ) psi

        """
        self.timer.start('Apply time-dependent operators')
        self.td_overlap.update_k_point_projections(self.wfs, self.kpt, psi)
        self.td_hamiltonian.apply(self.kpt, psi, self.hpsit,
                                  calculate_P_ani=False)
        self.td_hamiltonian.calculate_paw_correction(psi, self.hpsit, self.wfs, self.kpt, self.v_at, calculate_P_ani=False)
        self.td_overlap.apply(psi, self.spsit, self.wfs, self.kpt, calculate_P_ani=False)
        self.timer.stop('Apply time-dependent operators')

        # psin[:] = self.spsit + .5J * self.time_step * self.hpsit
        psin[:] = self.spsit
        self.mblas.multi_zaxpy(.5j*self.time_step, self.hpsit, psin, len(psi))
        # Apply shift -i eps S t/2 
        #self.mblas.multi_zaxpy(.5j * self.time_step * (-self.shift), self.spsit, psin, len(psi))

class EhrenfestHGHSICN(ExplicitCrankNicolson):
    """Semi-implicit Crank-Nicolson propagator for Ehrenfest dynamics
       using HGH pseudopotentials
    
    """
    def __init__(self, td_density, td_hamiltonian, td_overlap, 
                 solver, preconditioner, gd, timer):
        """Create SemiImplicitCrankNicolson-object.
        
        Parameters
        ----------
        td_density: TimeDependentDensity
            the time-dependent density
        td_hamiltonian: TimeDependentHamiltonian
            the time-dependent hamiltonian
        td_overlap: TimeDependentOverlap
            the time-dependent overlap operator
        solver: LinearSolver
            solver for linear equations
        preconditioner: Preconditioner
            preconditioner for linear equations
        gd: GridDescriptor
            coarse (wavefunction) grid descriptor
        timer: Timer
            timer
        
        """
        ExplicitCrankNicolson.__init__(self, td_density, td_hamiltonian,
                          td_overlap, solver, preconditioner, gd, timer)

        self.old_kpt_u = None

    def propagate(self, time, time_step):
        """Propagate wavefunctions once.
        
        Parameters
        ----------
        time: float
            the current time
        time_step: float
            time step
        """

        self.niter = 0

        # Allocate old/temporary wavefunctions
        if self.old_kpt_u is None:
            self.old_kpt_u = []
            for kpt in self.wfs.kpt_u:
                old_kpt = DummyKPoint()
                old_kpt.psit_nG = self.gd.empty(n=len(kpt.psit_nG),
                                                dtype=complex)
                self.old_kpt_u.append(old_kpt)

        if self.tmp_kpt_u is None:
            self.tmp_kpt_u = []
            for kpt in self.wfs.kpt_u:
                tmp_kpt = DummyKPoint()
                tmp_kpt.psit_nG = self.gd.empty(n=len(kpt.psit_nG),
                                                dtype=complex)
                self.tmp_kpt_u.append(tmp_kpt)

        # Allocate memory for Crank-Nicolson stuff
        nvec = len(self.wfs.kpt_u[0].psit_nG)
        if self.hpsit is None:
            self.hpsit = self.gd.zeros(nvec, dtype=complex)
        if self.spsit is None:
            self.spsit = self.gd.zeros(nvec, dtype=complex)


        self.timer.start('Update time-dependent operators')

        # Update overlap S(t) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t) based on the wavefunctions psit_nG
        # in kpt_u for t = time. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Update Hamiltonian H(t) to reflect density rho(t)
        self.td_hamiltonian.update(self.td_density.get_density(), time)

        self.timer.stop('Update time-dependent operators')

        # Copy current wavefunctions psit_nG to work and old wavefunction arrays
        for u, kpt in enumerate(self.wfs.kpt_u):
            self.old_kpt_u[u].psit_nG[:] = kpt.psit_nG
            self.tmp_kpt_u[u].psit_nG[:] = kpt.psit_nG

        #print 'P_ani[0] =', self.wfs.kpt_u[0].P_ani[0]
        #print self.test


        # Predictor step
        # Overwrite psit_nG in tmp_kpt_u by (1 - i S^(-1)(t) H(t) dt) psit_nG
        # from corresponding kpt_u in a Euler step before predicting psit(t+dt)
        for [kpt, rhs_kpt] in zip(self.wfs.kpt_u, self.tmp_kpt_u):
            self.solve_propagation_equation(kpt, rhs_kpt, time_step, guess=False)

        self.timer.start('Update time-dependent operators')

        # Update overlap S(t+dt) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t+dt) based on the wavefunctions psit_nG in
        # kpt_u for t = time+time_step. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Estimate Hamiltonian H(t+dt/2) by averaging H(t) and H(t+dt)
        # and retain the difference for a half-way Hamiltonian dH(t+dt/2).
        self.td_hamiltonian.half_update(self.td_density.get_density(),
                                        time + time_step)

        # Estimate overlap S(t+dt/2) by averaging S(t) and S(t+dt) #TODO!!!
        self.td_overlap.half_update(self.wfs)

        self.timer.stop('Update time-dependent operators')

        # Corrector step
        # Use predicted psit_nG in kpt_u as an initial guess, whereas the old 
        # wavefunction in old_kpt_u are used to calculate rhs based on psit(t)
        for [kpt, rhs_kpt] in zip(self.wfs.kpt_u, self.old_kpt_u):

            self.solve_propagation_equation(kpt, rhs_kpt, time_step)
 
        # update projections before exiting
        self.td_overlap.update(self.wfs)
       
        return self.niter


    # ( S + i H dt/2 ) psit(t+dt) = ( S - i H dt/2 ) psit(t)
    def solve_propagation_equation(self, kpt, rhs_kpt, time_step, calculate_P_ani = False, guess=False):

        # kpt is guess, rhs_kpt is used to calculate rhs and is overwritten
        nvec = len(rhs_kpt.psit_nG)

        assert kpt != rhs_kpt, 'Data race condition detected'
        assert len(kpt.psit_nG) == nvec, 'Incompatible lhs/rhs vectors'

        self.timer.start('Apply time-dependent operators')
        # Store H psi(t) as hpsit and S psit(t) as spsit
        self.td_overlap.update_k_point_projections(self.wfs, kpt, rhs_kpt.psit_nG)
        self.td_hamiltonian.apply(kpt, rhs_kpt.psit_nG, self.hpsit,
                                  calculate_P_ani=False)          
        self.td_overlap.apply(rhs_kpt.psit_nG, self.spsit, self.wfs, kpt,
                              calculate_P_ani=False)
        self.timer.stop('Apply time-dependent operators')

        # Update rhs_kpt.psit_nG to reflect ( S - i H dt/2 ) psit(t)
        #rhs_kpt.psit_nG[:] = self.spsit - .5J * self.hpsit * time_step
        rhs_kpt.psit_nG[:] = self.spsit
        self.mblas.multi_zaxpy(-.5j*time_step, self.hpsit, rhs_kpt.psit_nG, nvec)
        # Apply shift -i eps S t/2
        #self.mblas.multi_zaxpy(-.5j*time_step * (-self.shift), self.spsit, rhs_kpt.psit_nG, nvec)

        # Information needed by solver.solve -> self.dot
        self.kpt = kpt
        self.time_step = time_step

        # Solve A x = b where A is (S + i H dt/2) and b = rhs_kpt.psit_nG
        self.niter += self.solver.solve(self, kpt.psit_nG, rhs_kpt.psit_nG)

        # Apply shift exp(i eps t)
        #self.phase_shift = np.exp(1.0J * self.shift * time_step)
        #self.mblas.multi_scale(self.phase_shift, kpt.psit_nG, nvec)

    # ( S + i H dt/2 ) psi
    def dot(self, psi, psin):
        """Applies the propagator matrix to the given wavefunctions.

        Parameters
        ----------
        psi: List of coarse grids
            the known wavefunctions
        psin: List of coarse grids
            the result ( S + i H dt/2 ) psi

        """
        ExplicitCrankNicolson.dot(self, psi, psin)





###############################################################################
# EnforcedTimeReversalSymmetryCrankNicolson
###############################################################################
class EnforcedTimeReversalSymmetryCrankNicolson(ExplicitCrankNicolson):
    """Enforced time-reversal symmetry Crank-Nicolson propagator
    
    Crank-Nicolson propagator, which first approximates the time-dependent
    Hamiltonian to be unchanged during one iteration step to predict future
    wavefunctions. Then the approximations for the future wavefunctions are
    used to approximate the Hamiltonian in the future.
    
    (S(t) + .5j dt H(t) / hbar) psi(t+dt) = (S(t) - .5j dt H(t) / hbar) psi(t)
    (S(t) + .5j dt H(t+dt) / hbar) psi(t+dt) 
    = (S(t) - .5j dt H(t) / hbar) psi(t)
    
    """
    
    def __init__(self, td_density, td_hamiltonian, td_overlap, 
                 solver, preconditioner, gd, timer):
        """Create SemiImplicitCrankNicolson-object.
        
        Parameters
        ----------
        td_density: TimeDependentDensity
            the time-dependent density
        td_hamiltonian: TimeDependentHamiltonian
            the time-dependent hamiltonian
        td_overlap: TimeDependentOverlap
            the time-dependent overlap operator
        solver: LinearSolver
            solver for linear equations
        preconditioner: Preconditioner
            preconditioner for linear equations
        gd: GridDescriptor
            coarse (wavefunction) grid descriptor
        timer: Timer
            timer
        
        """
        ExplicitCrankNicolson.__init__(self, td_density, td_hamiltonian,
                          td_overlap, solver, preconditioner, gd, timer)

        self.old_kpt_u = None


    def propagate(self, time, time_step, update_callback=None):
        """Propagate wavefunctions once.
        
        Parameters
        ----------
        time: float
            the current time
        time_step: float
            time step
        """

        self.niter = 0

        # Allocate old/temporary wavefunctions
        if self.old_kpt_u is None:
            self.old_kpt_u = []
            for kpt in self.wfs.kpt_u:
                old_kpt = DummyKPoint()
                old_kpt.psit_nG = self.gd.empty(n=len(kpt.psit_nG),
                                                dtype=complex)
                self.old_kpt_u.append(old_kpt)

        if self.tmp_kpt_u is None:
            self.tmp_kpt_u = []
            for kpt in self.wfs.kpt_u:
                tmp_kpt = DummyKPoint()
                tmp_kpt.psit_nG = self.gd.empty(n=len(kpt.psit_nG),
                                                dtype=complex)
                self.tmp_kpt_u.append(tmp_kpt)

        # Allocate memory for Crank-Nicolson stuff
        nvec = len(self.wfs.kpt_u[0].psit_nG)
        if self.hpsit is None:
            self.hpsit = self.gd.zeros(nvec, dtype=complex)
        if self.spsit is None:
            self.spsit = self.gd.zeros(nvec, dtype=complex)


        self.timer.start('Update time-dependent operators')

        # Update overlap S(t) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t) based on the wavefunctions psit_nG
        # in kpt_u for t = time. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Update Hamiltonian H(t) to reflect density rho(t)
        self.td_hamiltonian.update(self.td_density.get_density(), time)

        self.timer.stop('Update time-dependent operators')

        # Copy current wavefunctions psit_nG to work and old wavefunction arrays
        for u, kpt in enumerate(self.wfs.kpt_u):
            self.old_kpt_u[u].psit_nG[:] = kpt.psit_nG
            self.tmp_kpt_u[u].psit_nG[:] = kpt.psit_nG


        # Predictor step
        for [kpt, rhs_kpt] in zip(self.wfs.kpt_u, self.old_kpt_u):
            self.create_rhs(rhs_kpt, kpt, time_step)

        if update_callback is not None:
            update_callback()

            
        for [kpt, rhs_kpt] in zip(self.wfs.kpt_u, self.old_kpt_u):            
            self.solve_propagation_equation(kpt, rhs_kpt, time_step)


        self.timer.start('Update time-dependent operators')

        # Update overlap S(t+dt) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t+dt) based on the wavefunctions psit_nG in
        # kpt_u for t = time+time_step. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Estimate Hamiltonian H(t+dt/2) by averaging H(t) and H(t+dt)
        # and retain the difference for a half-way Hamiltonian dH(t+dt/2).
        self.td_hamiltonian.update(self.td_density.get_density(),
                                   time + time_step)

        # Estimate overlap S(t+dt/2) by averaging S(t) and S(t+dt) #TODO!!!
        self.td_overlap.update(self.wfs)

        self.timer.stop('Update time-dependent operators')

        # Corrector step
        # Use predicted psit_nG in kpt_u as an initial guess, whereas the old 
        # wavefunction in old_kpt_u are used to calculate rhs based on psit(t)
        for [kpt, rhs_kpt] in zip(self.wfs.kpt_u, self.old_kpt_u):
            self.solve_propagation_equation(kpt, rhs_kpt, time_step)

        # update projections before exiting
        self.td_overlap.update(self.wfs)
        
        return self.niter

    # Create RHS
    def create_rhs(self, rhs_kpt, kpt, time_step):
        # kpt is guess, rhs_kpt is used to calculate rhs and is overwritten
        nvec = len(rhs_kpt.psit_nG)

        assert kpt != rhs_kpt, 'Data race condition detected'
        assert len(kpt.psit_nG) == nvec, 'Incompatible lhs/rhs vectors'

        self.timer.start('Apply time-dependent operators')
        # Store H psi(t) as hpsit and S psit(t) as spsit
        self.td_overlap.update_k_point_projections(self.wfs, kpt, rhs_kpt.psit_nG)
        self.td_hamiltonian.apply(kpt, rhs_kpt.psit_nG, self.hpsit,
                                  calculate_P_ani=False)
        self.td_overlap.apply(rhs_kpt.psit_nG, self.spsit, self.wfs, kpt,
                              calculate_P_ani=False)
        self.timer.stop('Apply time-dependent operators')

        # Update rhs_kpt.psit_nG to reflect ( S - i H dt/2 ) psit(t)
        #rhs_kpt.psit_nG[:] = self.spsit - .5J * self.hpsit * time_step
        rhs_kpt.psit_nG[:] = self.spsit
        self.mblas.multi_zaxpy(-.5j*time_step, self.hpsit, rhs_kpt.psit_nG, nvec)

    # ( S + i H(t+dt) dt/2 ) psit(t+dt) = ( S - i H(t) dt/2 ) psit(t)
    # rhs_kpt = ( S - i H(t) dt/2 ) psit(t)
    def solve_propagation_equation(self, kpt, rhs_kpt, time_step, guess=False):
        # Information needed by solver.solve -> self.dot
        self.kpt = kpt
        self.time_step = time_step

        # Solve A x = b where A is (S + i H dt/2) and b = rhs_kpt.psit_nG
        self.niter += self.solver.solve(self, kpt.psit_nG, rhs_kpt.psit_nG)

        # Apply shift exp(i eps t)
        #self.phase_shift = np.exp(1.0J * self.shift * time_step)
        #self.mblas.multi_scale(self.phase_shift, kpt.psit_nG, nvec)

    # ( S + i H dt/2 ) psi
    def dot(self, psi, psin):
        """Applies the propagator matrix to the given wavefunctions.

        Parameters
        ----------
        psi: List of coarse grids
            the known wavefunctions
        psin: List of coarse grids
            the result ( S + i H dt/2 ) psi

        """
        ExplicitCrankNicolson.dot(self, psi, psin)
        # Apply shift -i eps S t/2 
        #self.mblas.multi_zaxpy(.5j * self.time_step * (-self.shift), self.spsit, psin, len(psi))






###############################################################################
# AbsorptionKick
###############################################################################
class AbsorptionKick(ExplicitCrankNicolson):
    """Absorption kick propagator
    
    Absorption kick propagator::

      (S(t) + .5j dt p.r / hbar) psi(0+) = (S(t) - .5j dt p.r / hbar) psi(0-)

    where ``|p| = (eps e / hbar)``, and eps is field strength, e is elementary
    charge.
    
    """
    
    def __init__(self, wfs, abs_kick_hamiltonian, td_overlap,
                 solver, preconditioner, gd, timer):
        """Create AbsorptionKick-object.
        
        Parameters
        ----------
        wfs: FDWaveFunctions
            time-independent grid-based wavefunctions
        abs_kick_hamiltonian: AbsorptionKickHamiltonian
            the absorption kick hamiltonian
        td_overlap: TimeDependentOverlap
            the time-dependent overlap operator
        solver: LinearSolver
            solver for linear equations
        preconditioner: Preconditioner
            preconditioner for linear equations
        gd: GridDescriptor
            coarse (wavefunction) grid descriptor
        timer: Timer
            timer

        """
        ExplicitCrankNicolson.__init__(self, DummyDensity(wfs),
                        abs_kick_hamiltonian, td_overlap, solver,
                        preconditioner, gd, timer)


    def kick(self):
        """Excite all possible frequencies.
        
        """

        # if rank == 0:
        #     self.text('Kick iterations = ', self.td_hamiltonian.iterations)

        for l in range(self.td_hamiltonian.iterations):
            self.propagate(0, 1.0)
            # if rank == 0:
            #     self.text('.')
        # if rank == 0:
        #     print ''
        





###############################################################################
# SemiImpicitTaylorExponential
###############################################################################
class SemiImplicitTaylorExponential(DummyPropagator):
    """Semi-implicit Taylor exponential propagator 
    exp(-i S^-1 H t) = 1 - i S^-1 H t + (1/2) (-i S^-1 H t)^2 + ...
    
    """
    
    def __init__(self, td_density, td_hamiltonian, td_overlap, solver,
                 preconditioner, gd, timer, degree=4):
        """Create SemiImplicitTaylorExponential-object.
        
        Parameters
        ----------
        td_density: TimeDependentDensity
            the time-dependent density
        td_hamiltonian: TimeDependentHamiltonian
            the time-dependent hamiltonian
        td_overlap: TimeDependentOverlap
            the time-dependent overlap operator
        solver: LinearSolver
            solver for linear equations
        preconditioner: Preconditioner
            preconditioner
        gd: GridDescriptor
            coarse (wavefunction) grid descriptor
        timer: Timer
            timer
        degree: integer
            Degree of the Taylor polynomial (default is 4)
        
        """
        DummyPropagator.__init__(self, td_density, td_hamiltonian, td_overlap,
                                 solver, preconditioner, gd, timer)

        self.degree = degree

        self.tmp_kpt_u = None
        self.psin = None
        self.hpsit = None
        
        
    def propagate(self, time, time_step):
        """Propagate wavefunctions once.
        
        Parameters
        ----------
        time: float
            the current time
        time_step: float
            time step
        """

        self.niter = 0

        # Allocate temporary wavefunctions
        if self.tmp_kpt_u is None:
            self.tmp_kpt_u = []
            for kpt in self.wfs.kpt_u:
                tmp_kpt = DummyKPoint()
                tmp_kpt.psit_nG = self.gd.empty(n=len(kpt.psit_nG),
                                                dtype=complex)
                self.tmp_kpt_u.append(tmp_kpt)

        # Allocate memory for Taylor exponential stuff
        nvec = len(self.wfs.kpt_u[0].psit_nG)
        if self.psin is None:
            self.psin = self.gd.zeros(nvec, dtype=complex)
        if self.hpsit is None:
            self.hpsit = self.gd.zeros(nvec, dtype=complex)
        

        self.timer.start('Update time-dependent operators')

        # Update overlap S(t) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t) based on the wavefunctions psit_nG
        # in kpt_u for t = time. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Update Hamiltonian H(t) to reflect density rho(t)
        self.td_hamiltonian.update(self.td_density.get_density(), time)

        self.timer.stop('Update time-dependent operators')


        # copy current wavefunctions to temporary variable
        for u, kpt in enumerate(self.wfs.kpt_u):
            self.tmp_kpt_u[u].psit_nG[:] = kpt.psit_nG

        # predict for each k-point
        for kpt in self.wfs.kpt_u:
            self.solve_propagation_equation(kpt, time_step)


        self.timer.start('Update time-dependent operators')

        # Update overlap S(t+dt) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t+dt) based on the wavefunctions psit_nG in
        # kpt_u for t = time+time_step. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Estimate Hamiltonian H(t+dt/2) by averaging H(t) and H(t+dt)
        self.td_hamiltonian.half_update(self.td_density.get_density(),
                                        time + time_step)

        # Estimate overlap S(t+dt/2) by averaging S(t) and S(t+dt) #TODO!!!
        self.td_overlap.half_update(self.wfs)

        self.timer.stop('Update time-dependent operators')


        # propagate psit(t), not psit(t+dt), in correct
        for u, kpt in enumerate(self.wfs.kpt_u):
            kpt.psit_nG[:] = self.tmp_kpt_u[u].psit_nG

        # correct for each k-point
        for kpt in self.wfs.kpt_u:
            self.solve_propagation_equation(kpt, time_step)
        
        # update projections before exiting
        self.td_overlap.update(self.wfs)

        return self.niter

    # psi(t) = exp(-i t S^-1 H) psi(0)
    # psi(t) = 1  + (-i S^-1 H t) (1 + (1/2) (-i S^-1 H t) (1 + ... ) )
    def solve_propagation_equation(self, kpt, time_step):

        nvec = len(kpt.psit_nG)

        # Information needed by solver.solve -> self.dot
        self.kpt = kpt
        self.time_step = time_step
            
        # psin = psi(0)
        self.psin[:] = kpt.psit_nG
        for k in range(self.degree,0,-1):
            # psin = psi(0) + (1/k) (-i S^-1 H t) psin
            self.td_hamiltonian.apply(kpt, self.psin, self.hpsit)
            # S psin = H psin
            self.psin[:] = self.hpsit
            self.niter += self.solver.solve(self, self.psin, self.hpsit)
            #print 'Linear solver iterations = ', self.solver.iterations
            # psin = psi(0) + (-it/k) S^-1 H psin
            self.mblas.multi_scale(-1.0j*time_step/k, self.psin, nvec)
            self.mblas.multi_zaxpy(1.0, kpt.psit_nG, self.psin, nvec)

        kpt.psit_nG[:] = self.psin

    def dot(self, psit, spsit):
        self.td_overlap.apply(psit, spsit, self.wfs, self.kpt)




###############################################################################
# SemiImplicitKrylovExponential
###############################################################################
class SemiImplicitKrylovExponential(DummyPropagator):
    """Semi-implicit Krylov exponential propagator
    
    
    """
    
    def __init__(self, td_density, td_hamiltonian, td_overlap, solver,
                 preconditioner, gd, timer, degree=4):
        """Create SemiImplicitKrylovExponential-object.
        
        Parameters
        ----------
        td_density: TimeDependentDensity
            the time-dependent density
        td_hamiltonian: TimeDependentHamiltonian
            the time-dependent hamiltonian
        td_overlap: TimeDependentOverlap
            the time-dependent overlap operator
        solver: LinearSolver
            solver for linear equations
        preconditioner: Preconditioner
            preconditioner
        gd: GridDescriptor
            coarse (wavefunction) grid descriptor
        timer: Timer
            timer
        degree: integer
            Degree of the Krylov subspace (default is 4)
        
        """
        DummyPropagator.__init__(self, td_density, td_hamiltonian, td_overlap,
                                 solver, preconditioner, gd, timer)

        self.kdim = degree + 1
        
        self.tmp_kpt_u = None
        self.lm = None
        self.em = None
        self.hm = None
        self.sm = None
        self.xm = None
        self.qm = None
        self.Hqm = None
        self.Sqm = None
        self.rqm = None
        
        
    def propagate(self, time, time_step):
        """Propagate wavefunctions once.
        
        Parameters
        ----------
        time: float
            the current time
        time_step: float
            time step
        """

        self.niter = 0

        # Allocate temporary wavefunctions
        if self.tmp_kpt_u is None:
            self.tmp_kpt_u = []
            for kpt in self.wfs.kpt_u:
                tmp_kpt = DummyKPoint()
                tmp_kpt.psit_nG = self.gd.empty(n=len(kpt.psit_nG),
                                                dtype=complex)
                self.tmp_kpt_u.append(tmp_kpt)


        # Allocate memory for Krylov subspace stuff
        nvec = len(self.wfs.kpt_u[0].psit_nG)

        # em = (wfs, degree)
        if self.em is None:
            self.em = np.zeros((nvec, self.kdim), float)

        # lm = (wfs)
        if self.lm is None:
            self.lm = np.zeros((nvec,), complex)

        # hm = (wfs, degree, degree)
        if self.hm is None:
            self.hm = np.zeros((nvec, self.kdim, self.kdim), complex)
        # sm = (wfs, degree, degree)
        if self.sm is None:
            self.sm = np.zeros((nvec, self.kdim, self.kdim), complex)
        # xm = (wfs, degree, degree)
        if self.xm is None:
            self.xm = np.zeros((nvec, self.kdim, self.kdim), complex)

        # qm = (degree, wfs, nx, ny, nz) 
        if self.qm is None:
            self.qm = self.gd.zeros((self.kdim, nvec), dtype=complex)
        # H qm = (degree, wfs, nx, ny, nz) 
        if self.Hqm is None:
            self.Hqm = self.gd.zeros((self.kdim, nvec), dtype=complex)
        # S qm = (degree, wfs, nx, ny, nz) 
        if self.Sqm is None:
            self.Sqm = self.gd.zeros((self.kdim, nvec), dtype=complex)
        # rqm = (wfs, nx, ny, nz) 
        if self.rqm is None:
            self.rqm = self.gd.zeros((nvec,), dtype=complex)

        
        self.time_step = time_step


        self.timer.start('Update time-dependent operators')

        # Update overlap S(t) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t) based on the wavefunctions psit_nG
        # in kpt_u for t = time. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Update Hamiltonian H(t) to reflect density rho(t)
        self.td_hamiltonian.update(self.td_density.get_density(), time)

        self.timer.stop('Update time-dependent operators')


        # copy current wavefunctions to temporary variable
        for u, kpt in enumerate(self.wfs.kpt_u):
            self.tmp_kpt_u[u].psit_nG[:] = kpt.psit_nG

        # predict for each k-point
        for kpt in self.wfs.kpt_u:
            self.solve_propagation_equation(kpt, time_step)


        self.timer.start('Update time-dependent operators')

        # Update overlap S(t+dt) of kpt.psit_nG in kpt.P_ani.
        self.td_overlap.update(self.wfs)

        # Calculate density rho(t+dt) based on the wavefunctions psit_nG in
        # kpt_u for t = time+time_step. Updates wfs.D_asp based on kpt.P_ani.
        self.td_density.update()

        # Estimate Hamiltonian H(t+dt/2) by averaging H(t) and H(t+dt)
        self.td_hamiltonian.half_update(self.td_density.get_density(),
                                        time + time_step)

        # Estimate overlap S(t+dt/2) by averaging S(t) and S(t+dt) #TODO!!!
        self.td_overlap.half_update(self.wfs)

        self.timer.stop('Update time-dependent operators')


        # propagate psit(t), not psit(t+dt), in correct
        for u, kpt in enumerate(self.wfs.kpt_u):
            kpt.psit_nG[:] = self.tmp_kpt_u[u].psit_nG

        # correct for each k-point
        for kpt in self.wfs.kpt_u:
            self.solve_propagation_equation(kpt, time_step)
        
        # update projections before exiting
        self.td_overlap.update(self.wfs)

        return self.niter
    
    # psi(t) = exp(-i t S^-1 H) psi(0)
    def solve_propagation_equation(self, kpt, time_step):

        nvec = len(kpt.psit_nG)
        tmp = np.zeros((nvec,), complex)
        xm_tmp = np.zeros((nvec, self.kdim), complex)

        qm = self.qm
        Hqm = self.Hqm
        Sqm = self.Sqm

        # Information needed by solver.solve -> self.dot
        self.kpt = kpt

        scale = self.create_krylov_subspace(kpt,
                                            self.td_hamiltonian,
                                            self.td_overlap,
                                            self.qm, self.Hqm, self.Sqm)

        # Calculate hm and sm
        for i in range(self.kdim):
            for j in range(self.kdim):
                self.mblas.multi_zdotc(tmp, qm[i], Hqm[j], nvec)
                tmp *= self.gd.dv
                for k in range(nvec):
                    self.hm[k][i][j] = tmp[k]
                self.mblas.multi_zdotc(tmp, qm[i], Sqm[j], nvec)
                tmp *= self.gd.dv
                for k in range(nvec):
                    self.sm[k][i][j] = tmp[k]

        #print 'Hm ='
        #print np.round(self.hm*1e4) / 1e4
        #print 'log Hm ='
        #print np.round(np.log(self.hm)*1e2)/1e2
        #print 'Sm ='
        #print np.round(self.sm*1e4) / 1e4

        #print np.round(np.log10(np.abs(np.linalg.eigh(self.hm[0])[1]))*1e6)/1e6

        # Diagonalize
        # Propagate
        # psi(t) = Qm Xm exp(-i Em t) Xm^H Sm e_1
        #        = Qm Xm exp(-i Em t) Sm Qm^H S psi(0) ???
        #        = Qm Xm exp(-i Em t) y
        #        = Qm Xm z
        # y = Sm Qm^H S psi(0) = Xm^H Sm e_1
        # if Sm = I then y is the first row of Xm^*
        # and z = exp(-i Em t) y
        for k in range(nvec):
            (self.em[k], self.xm[k]) = np.linalg.eigh(self.hm[k])
        #print 'Em = ', self.em
        #for k in range(nvec):
            #print 'Xm',k,' = '
            #print self.xm[k]

        #print self.em[0] * (-1.0J*self.time_step)
        self.em = np.exp(self.em * (-1.0j*time_step))
        #print self.em[0]
        #print np.linalg.eigh(self.hm[0])
        for k in range(nvec):
            z = self.em[k] * np.conj(self.xm[k,0])
            xm_tmp[k][:] = np.dot(self.xm[k], z)
        #print xm_tmp
        kpt.psit_nG[:] = 0.0
        for k in range(nvec):
            for i in range(self.kdim):
                #print 'Xm_tmp[',k,'][',i,'] = ', xm_tmp[k][i]
                axpy(xm_tmp[k][i] / scale[k], 
                     self.qm[i][k], kpt.psit_nG[k])

        #print self.qm
        #print kpt.psit_nG


    # Create Krylov subspace
    #    K_v = { psi, S^-1 H psi, (S^-1 H)^2 psi, ... }
    def create_krylov_subspace(self, kpt, h, s, qm, Hqm, Sqm):
        nvec = len(kpt.psit_nG)
        # tmp = (wfs)
        tmp = np.zeros((nvec,), complex)
        scale = np.zeros((nvec,), complex)
        scale[:] = 0.0
        rqm = self.rqm

        # q_0 = psi
        rqm[:] = kpt.psit_nG

        for i in range(self.kdim):
            qm[i][:] = rqm

            # S orthogonalize
            # q_i = q_i - sum_j<i <q_j|S|q_i> q_j
            for j in range(i):
                self.mblas.multi_zdotc(tmp, qm[i], Sqm[j], nvec)
                tmp *= self.gd.dv
                tmp = np.conj(tmp)
                self.mblas.multi_zaxpy(-tmp, qm[j], qm[i], nvec)

            # S q_i
            s.apply(qm[i], Sqm[i], self.wfs, kpt)
            self.mblas.multi_zdotc(tmp, qm[i], Sqm[i], nvec)
            tmp *= self.gd.dv
            self.mblas.multi_scale(1./np.sqrt(tmp), qm[i], nvec)
            self.mblas.multi_scale(1./np.sqrt(tmp), Sqm[i], nvec)
            if i == 0:
                scale[:] = 1/np.sqrt(tmp)
                #print 'Scale', scale

            # H q_i
            h.apply(kpt, qm[i], Hqm[i])

            # S r = H q_i, (if stuff, to save one inversion)
            if i+1 < self.kdim:
                rqm[:] = Hqm[i]
                self.solver.solve(self, rqm, Hqm[i])
                #print 'Linear solver iterations = ', self.solver.iterations

        #print '---'
        return scale

    def dot(self, psit, spsit):
        self.td_overlap.apply(psit, spsit, self.wfs, self.kpt)

    ### Below this, just for testing & debug
    def Sdot(self, psit, spsit):
        self.apply_preconditioner(psit, self.tmp)
        self.td_overlap.apply(self.tmp, spsit, self.wfs, self.kpt)

    def Hdot(self, psit, spsit):
        self.apply_preconditioner(psit, self.tmp)
        self.td_hamiltonian.apply(self.kpt, self.tmp, spsit)

    def inverse_overlap(self, kpt_u, degree):
        self.dot = self.Sdot
        self.kpt = kpt_u[0]
        nvec = len(self.kpt.psit_nG)
        nrm2 = np.zeros(nvec, dtype=complex)
        self.tmp = self.gd.zeros(n=nvec, dtype=complex)

        for i in range(10):
            self.solver.solve(self, self.kpt.psit_nG, self.kpt.psit_nG)
            self.mblas.multi_zdotc(nrm2, self.kpt.psit_nG, self.kpt.psit_nG, nvec)
            nrm2 *= self.gd.dv
            self.mblas.multi_scale(1/np.sqrt(nrm2), self.kpt.psit_nG, nvec)
        self.td_overlap.apply(self.kpt.psit_nG, self.tmp, self.wfs, self.kpt)
        self.mblas.multi_zdotc(nrm2, self.kpt.psit_nG, self.tmp, nvec)
        nrm2 *= self.gd.dv
        print('S min eig = ', nrm2)


    def overlap(self, kpt_u, degree):
        self.dot = self.Sdot
        self.kpt = kpt_u[0]
        nvec = len(self.kpt.psit_nG)
        nrm2 = np.zeros(nvec, dtype=complex)
        self.tmp = self.gd.zeros(n=nvec, dtype=complex)

        for i in range(100):
            self.tmp[:] = self.kpt.psit_nG
            self.td_overlap.apply(self.tmp, self.kpt.psit_nG, self.wfs, self.kpt)
            self.mblas.multi_zdotc(nrm2, self.kpt.psit_nG, self.kpt.psit_nG, nvec)
            nrm2 *= self.gd.dv
            self.mblas.multi_scale(1/np.sqrt(nrm2), self.kpt.psit_nG, nvec)
        self.td_overlap.apply(self.kpt.psit_nG, self.tmp, self.wfs, self.kpt)
        self.mblas.multi_zdotc(nrm2, self.kpt.psit_nG, self.tmp, nvec)
        nrm2 *= self.gd.dv
        print('S max eig = ', nrm2)


    def inverse_hamiltonian(self, kpt_u, degree):
        self.dot = self.Hdot
        self.kpt = kpt_u[0]
        nvec = len(self.kpt.psit_nG)
        nrm2 = np.zeros(nvec, dtype=complex)
        self.tmp = self.gd.zeros(n=nvec, dtype=complex)

        for i in range(10):
            self.solver.solve(self, self.kpt.psit_nG, self.kpt.psit_nG)
            self.mblas.multi_zdotc(nrm2, self.kpt.psit_nG, self.kpt.psit_nG, nvec)
            nrm2 *= self.gd.dv
            self.mblas.multi_scale(1/np.sqrt(nrm2), self.kpt.psit_nG, nvec)
        self.td_hamiltonian.apply(self.kpt, self.kpt.psit_nG, self.tmp)
        self.mblas.multi_zdotc(nrm2, self.kpt.psit_nG, self.tmp, nvec)
        nrm2 *= self.gd.dv
        print('H min eig = ', nrm2)


    def hamiltonian(self, kpt_u, degree):
        self.dot = self.Hdot
        self.kpt = kpt_u[0]
        nvec = len(self.kpt.psit_nG)
        nrm2 = np.zeros(nvec, dtype=complex)
        self.tmp = self.gd.zeros(n=nvec, dtype=complex)

        for i in range(100):
            self.tmp[:] = self.kpt.psit_nG
            self.td_hamiltonian.apply(self.kpt, self.tmp, self.kpt.psit_nG)
            self.mblas.multi_zdotc(nrm2, self.kpt.psit_nG, self.kpt.psit_nG, nvec)
            nrm2 *= self.gd.dv
            self.mblas.multi_scale(1/np.sqrt(nrm2), self.kpt.psit_nG, nvec)
        self.td_hamiltonian.apply(self.kpt, self.kpt.psit_nG, self.tmp)
        self.mblas.multi_zdotc(nrm2, self.kpt.psit_nG, self.tmp, nvec)
        nrm2 *= self.gd.dv
        print('H max eig = ', nrm2)
