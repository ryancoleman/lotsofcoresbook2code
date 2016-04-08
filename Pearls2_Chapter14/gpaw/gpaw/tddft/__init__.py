"""This module implements a class for (true) time-dependent density
functional theory calculations.

"""
from __future__ import print_function
import sys
import time
from math import log

import numpy as np
from ase.units import Bohr, Hartree

import gpaw.io
import gpaw.mpi as mpi
from gpaw.aseinterface import GPAW
from gpaw.mixer import DummyMixer
from gpaw.version import version
from gpaw.preconditioner import Preconditioner
from gpaw.lfc import LocalizedFunctionsCollection as LFC
from gpaw.tddft.units import attosec_to_autime, autime_to_attosec, \
                             eV_to_aufrequency, aufrequency_to_eV
from gpaw.tddft.utils import MultiBlas
from gpaw.tddft.bicgstab import BiCGStab
from gpaw.tddft.cscg import CSCG
from gpaw.tddft.propagators import \
    ExplicitCrankNicolson, \
    SemiImplicitCrankNicolson, \
    EhrenfestPAWSICN,\
    EhrenfestHGHSICN,\
    EnforcedTimeReversalSymmetryCrankNicolson, \
    SemiImplicitTaylorExponential, \
    SemiImplicitKrylovExponential, \
    AbsorptionKick
from gpaw.tddft.tdopers import \
    TimeDependentHamiltonian, \
    TimeDependentOverlap, \
    TimeDependentWaveFunctions, \
    TimeDependentDensity, \
    AbsorptionKickHamiltonian
from gpaw.tddft.abc import \
    LinearAbsorbingBoundary, \
    PML, \
    P4AbsorbingBoundary


# T^-1
# Bad preconditioner
class KineticEnergyPreconditioner:
    def __init__(self, gd, kin, dtype):
        self.preconditioner = Preconditioner(gd, kin, dtype)
        self.preconditioner.allocate()

    def apply(self, kpt, psi, psin):
        for i in range(len(psi)):
            psin[i][:] = self.preconditioner(psi[i], kpt.phase_cd, None, None)

# S^-1
class InverseOverlapPreconditioner:
    """Preconditioner for TDDFT."""
    def __init__(self, overlap):
        self.overlap = overlap

    def apply(self, kpt, psi, psin):
        self.overlap.apply_inverse(psi, psin, kpt)
# ^^^^^^^^^^


###########################
# Main class
###########################
class TDDFT(GPAW):
    """Time-dependent density functional theory calculation based on GPAW.
    
    This class is the core class of the time-dependent density functional
    theory implementation and is the only class which a user has to use.
    """
    
    def __init__(self, filename, td_potential=None, propagator='SICN', calculate_energy=True,
                 propagator_kwargs=None, solver='CSCG', tolerance=1e-8,
                 **kwargs):
        """Create TDDFT-object.
        
        Parameters:

        filename: string
            File containing ground state or time-dependent state to propagate
        td_potential: class, optional
            Function class for the time-dependent potential. Must have a method
            'strength(time)' which returns the strength of the linear potential
            to each direction as a vector of three floats.
        propagator:  {'SICN','ETRSCN','ECN','SITE','SIKE4','SIKE5','SIKE6'}
            Name of the time propagator for the Kohn-Sham wavefunctions
        solver: {'CSCG','BiCGStab'}
            Name of the iterative linear equations solver for time propagation
        tolerance: float
            Tolerance for the linear solver

        The following parameters can be used: `txt`, `parallel`, `communicator`
        `mixer` and `dtype`. The internal parameters `mixer` and `dtype` are
        strictly used to specify a dummy mixer and complex type respectively.
        """

        # Set initial time
        self.time = 0.0

        # Set initial kick strength
        self.kick_strength = np.array([0.0, 0.0, 0.0], dtype=float)

        # Set initial value of iteration counter
        self.niter = 0

        # Override default `mixer` and `dtype` given in InputParameters
        kwargs.setdefault('mixer', DummyMixer())
        kwargs.setdefault('dtype', complex)

        # Parallelization dictionary should also default to strided bands
        parallel = kwargs.setdefault('parallel', {})
        parallel.setdefault('stridebands', True)

        # Initialize paw-object without density mixing
        # NB: TDDFT restart files contain additional information which
        #     will override the initial settings for time/kick/niter.
        GPAW.__init__(self, filename, **kwargs)

        assert isinstance(self.wfs, TimeDependentWaveFunctions)
        assert isinstance(self.wfs.overlap, TimeDependentOverlap)

        # Prepare for dipole moment file handle
        self.dm_file = None

        # Initialize wavefunctions and density
        # (necessary after restarting from file)
        self.set_positions()

        # Don't be too strict
        self.density.charge_eps = 1e-5

        wfs = self.wfs
        self.rank = wfs.world.rank

        self.text('')
        self.text('')
        self.text('------------------------------------------')
        self.text('  Time-propagation TDDFT                  ')
        self.text('------------------------------------------')
        self.text('')

        self.text('Charge epsilon: ', self.density.charge_eps)

        # Time-dependent variables and operators
        self.td_potential = td_potential
        self.td_hamiltonian = TimeDependentHamiltonian(self.wfs, self.atoms,
                                  self.hamiltonian, td_potential)
        self.td_overlap = self.wfs.overlap #TODO remove this property
        self.td_density = TimeDependentDensity(self)

        # Solver for linear equations
        self.text('Solver: ', solver)
        if solver == 'BiCGStab':
            self.solver = BiCGStab(gd=wfs.gd, timer=self.timer,
                                   tolerance=tolerance)
        elif solver == 'CSCG':
            self.solver = CSCG(gd=wfs.gd, timer=self.timer,
                               tolerance=tolerance)
        else:
            raise RuntimeError('Solver %s not supported.' % solver)

        # Preconditioner
        # No preconditioner as none good found
        self.text('Preconditioner: ', 'None')
        self.preconditioner = None #TODO! check out SSOR preconditioning
        #self.preconditioner = InverseOverlapPreconditioner(self.overlap)
        #self.preconditioner = KineticEnergyPreconditioner(wfs.gd, self.td_hamiltonian.hamiltonian.kin, np.complex)

        # Time propagator
        self.text('Propagator: ', propagator)
        if propagator_kwargs is None:
            propagator_kwargs = {}
        if propagator == 'ECN':
            self.propagator = ExplicitCrankNicolson(self.td_density,
                self.td_hamiltonian, self.td_overlap, self.solver,
                self.preconditioner, wfs.gd, self.timer, **propagator_kwargs)
        elif propagator == 'SICN':
            self.propagator = SemiImplicitCrankNicolson(self.td_density,
                self.td_hamiltonian, self.td_overlap, self.solver,
                self.preconditioner, wfs.gd, self.timer, **propagator_kwargs)
        elif propagator == 'EFSICN':
            self.propagator = EhrenfestPAWSICN(self.td_density,
                self.td_hamiltonian, self.td_overlap, self.solver,
                self.preconditioner, wfs.gd, self.timer, **propagator_kwargs)
        elif propagator == 'EFSICN_HGH':
            self.propagator = EhrenfestHGHSICN(self.td_density,
                self.td_hamiltonian, self.td_overlap, self.solver,
                self.preconditioner, wfs.gd, self.timer, **propagator_kwargs)
        elif propagator == 'ETRSCN':
            self.propagator = EnforcedTimeReversalSymmetryCrankNicolson(
                self.td_density,
                self.td_hamiltonian, self.td_overlap, self.solver,
                self.preconditioner, wfs.gd, self.timer, **propagator_kwargs)
        elif propagator == 'SITE':
            self.propagator = SemiImplicitTaylorExponential(self.td_density,
                self.td_hamiltonian, self.td_overlap, self.solver,
                self.preconditioner, wfs.gd, self.timer, **propagator_kwargs)
        elif propagator == 'SIKE':
            self.propagator = SemiImplicitKrylovExponential(self.td_density,
                self.td_hamiltonian, self.td_overlap, self.solver,
                self.preconditioner, wfs.gd, self.timer, **propagator_kwargs)
        elif propagator.startswith('SITE') or propagator.startswith('SIKE'):
            raise DeprecationWarning('Use propagator_kwargs to specify degree.')
        else:
            raise RuntimeError('Time propagator %s not supported.' % propagator)

        if self.rank == 0:
            if wfs.kd.comm.size > 1:
                if wfs.nspins == 2:
                    self.text('Parallelization Over Spin')

                if wfs.gd.comm.size > 1:
                    self.text('Using Domain Decomposition: %d x %d x %d' %
                              tuple(wfs.gd.parsize_c))

                if wfs.band_comm.size > 1:
                    self.text('Parallelization Over bands on %d Processors' %
                              wfs.band_comm.size)
            self.text('States per processor = ', wfs.bd.mynbands)

        self.hpsit = None
        self.eps_tmp = None
        self.mblas = MultiBlas(wfs.gd)
        
        # Restarting an FDTD run generates hamiltonian.fdtd_poisson, which now overwrites hamiltonian.poisson
        if hasattr(self.hamiltonian, 'fdtd_poisson'):
            self.hamiltonian.poisson = self.hamiltonian.fdtd_poisson
            self.hamiltonian.poisson.set_grid_descriptor(self.density.finegd)

        # For electrodynamics mode
        if self.hamiltonian.poisson.get_description() == 'FDTD+TDDFT':
            self.initialize_FDTD()
            self.hamiltonian.poisson.print_messages(self.text)
            self.txt.flush()



        self.calculate_energy = calculate_energy
        if self.hamiltonian.xc.name.startswith('GLLB'):
             self.text("GLLB model potential. Not updating energy.")
             self.calculate_energy = False


    # Electrodynamics requires extra care
    def initialize_FDTD(self):
        
        # Sanity check
        assert(self.hamiltonian.poisson.get_description() == 'FDTD+TDDFT')
        
        self.hamiltonian.poisson.set_density(self.density)

        # The propagate calculation_mode causes classical part to evolve
        # in time when self.hamiltonian.poisson.solve(...) is called
        self.hamiltonian.poisson.set_calculation_mode('propagate')

        # During each time step, self.hamiltonian.poisson.solve may be called several
        # times (depending on the used propagator). Using the attached observer one
        # ensures that actual propagation takes place only once. This is because
        # the FDTDPoissonSolver changes the calculation_mode from propagate to
        # something else when the propagation is finished.
        self.attach(self.hamiltonian.poisson.set_calculation_mode, 1, 'propagate')


    def set(self, **kwargs):
        p = self.input_parameters

        # Special treatment for dictionary parameters:
        for name in ['parallel']:
            if kwargs.get(name) is not None:
                tmp = p[name]
                for key in kwargs[name]:
                    if not key in tmp:
                        raise KeyError('Unknown subparameter "%s" in '
                                       'dictionary parameter "%s"' % (key,
                                                                      name))
                tmp.update(kwargs[name])
                kwargs[name] = tmp

        for key in kwargs:
            # Only whitelisted arguments can be changed after initialization
            if self.initialized and key not in ['txt']:
                raise TypeError("Keyword argument '%s' is immutable." % key)

            if key in ['txt', 'parallel', 'communicator']:
                continue
            elif key == 'mixer':
                if not isinstance(kwargs[key], DummyMixer):
                    raise ValueError("Mixer must be of type DummyMixer.")
            elif key == 'dtype':
                if kwargs[key] is not complex:
                    raise ValueError("TDDFT calculation must be complex.")
            elif key in ['parsize', 'parsize_bands', 'parstride_bands']:
                name = {'parsize': 'domain',
                        'parsize_bands': 'band',
                        'parstride_bands': 'stridebands'}[key]
                raise DeprecationWarning(
                    'Keyword argument has been moved ' +
                    "to the 'parallel' dictionary keyword under '%s'." % name)
            else:
                raise TypeError("Unknown keyword argument: '%s'" % key)

        p.update(kwargs)

    def read(self, reader, read_projections=True):
        assert reader.has_array('PseudoWaveFunctions')
        GPAW.read(self, reader, read_projections)

    def propagate(self, time_step, iterations, dipole_moment_file=None,
                  restart_file=None, dump_interval=100):
        """Propagates wavefunctions.
        
        Parameters:

        time_step: float
            Time step in attoseconds (10^-18 s), e.g., 4.0 or 8.0
        iterations: integer
            Iterations, e.g., 20 000 as / 4.0 as = 5000
        dipole_moment_file: string, optional
            Name of the data file where to the time-dependent dipole
            moment is saved
        restart_file: string, optional
            Name of the restart file
        dump_interval: integer
            After how many iterations restart data is dumped
        
        """

        if self.rank == 0:
            self.text()
            self.text('Starting time: %7.2f as'
                      % (self.time * autime_to_attosec))
            self.text('Time step:     %7.2f as' % time_step)
            header = """\
                        Simulation      Total        log10     Iterations:
             Time          time         Energy       Norm      Propagator"""
            self.text()
            self.text(header)


        # Convert to atomic units
        time_step = time_step * attosec_to_autime
        
        if dipole_moment_file is not None:
            self.initialize_dipole_moment_file(dipole_moment_file)

        # Set these as class properties for use of observers
        self.time_step = time_step
        self.dump_interval = dump_interval

        niterpropagator = 0
        maxiter = self.niter + iterations

        # Let FDTD part know the time step
        if self.hamiltonian.poisson.get_description() == 'FDTD+TDDFT':
            self.hamiltonian.poisson.set_time(self.time)
            self.hamiltonian.poisson.set_time_step(self.time_step)
            
        self.timer.start('Propagate')
        while self.niter < maxiter:
            norm = self.density.finegd.integrate(self.density.rhot_g)

            # Write dipole moment at every iteration
            if dipole_moment_file is not None:
                self.update_dipole_moment_file(norm)

            # print output (energy etc.) every 10th iteration
            if self.niter % 10 == 0:
                self.get_td_energy()
                
                T = time.localtime()
                if self.rank == 0:
                    iter_text = 'iter: %3d  %02d:%02d:%02d %11.2f' \
                                '   %13.6f %9.1f %10d'
                    self.text(iter_text %
                              (self.niter, T[3], T[4], T[5],
                               self.time * autime_to_attosec,
                               self.Etot, log(abs(norm)+1e-16)/log(10),
                               niterpropagator))

                    self.txt.flush()


            # Propagate the Kohn-Shame wavefunctions a single timestep
            niterpropagator = self.propagator.propagate(self.time, time_step)
            self.time += time_step
            self.niter += 1

            # Call registered callback functions
            self.call_observers(self.niter)

            # Write restart data
            if restart_file is not None and self.niter % dump_interval == 0:
                self.write(restart_file, 'all')
                if self.rank == 0:
                    print('Wrote restart file.')
                    print(self.niter, ' iterations done. Current time is ', \
                        self.time * autime_to_attosec, ' as.')

        self.timer.stop('Propagate')

        # Write final results and close dipole moment file
        if dipole_moment_file is not None:
            #TODO final iteration is propagated, but nothing is updated
            #norm = self.density.finegd.integrate(self.density.rhot_g)
            #self.finalize_dipole_moment_file(norm)
            self.finalize_dipole_moment_file()

        # Finalize FDTDPoissonSolver
        if self.hamiltonian.poisson.get_description() == 'FDTD+TDDFT':
            self.hamiltonian.poisson.finalize_propagation()

        # Call registered callback functions
        self.call_observers(self.niter, final=True)

        if restart_file is not None:
            self.write(restart_file, 'all')

    def initialize_dipole_moment_file(self, dipole_moment_file):
        if self.rank == 0:
            if self.dm_file is not None and not self.dm_file.closed:
                raise RuntimeError('Dipole moment file is already open')

            if self.time == 0.0:
                mode = 'w'
            else:
                # We probably continue from restart
                mode = 'a'

            self.dm_file = file(dipole_moment_file, mode)

            # If the dipole moment file is empty, add a header
            if self.dm_file.tell() == 0:
                header = '# Kick = [%22.12le, %22.12le, %22.12le]\n' \
                    % (self.kick_strength[0], self.kick_strength[1], \
                       self.kick_strength[2])
                header += '# %15s %15s %22s %22s %22s\n' \
                    % ('time', 'norm', 'dmx', 'dmy', 'dmz')
                self.dm_file.write(header)
                self.dm_file.flush()
        

    def update_dipole_moment_file(self, norm):
        dm = self.density.finegd.calculate_dipole_moment(self.density.rhot_g)
        
        if self.hamiltonian.poisson.get_description() == 'FDTD+TDDFT':
            dm += self.hamiltonian.poisson.get_classical_dipole_moment()

        if self.rank == 0:
            line = '%20.8lf %20.8le %22.12le %22.12le %22.12le\n' \
                % (self.time, norm, dm[0], dm[1], dm[2])
            self.dm_file.write(line)
            self.dm_file.flush()

    def finalize_dipole_moment_file(self, norm=None):
        if norm is not None:
            self.update_dipole_moment_file(norm)

        if self.rank == 0:
            self.dm_file.close()
            self.dm_file = None

    def update_eigenvalues(self):

        kpt_u = self.wfs.kpt_u
        if self.hpsit is None:
            self.hpsit = self.wfs.gd.zeros(len(kpt_u[0].psit_nG),
                                           dtype=complex)
        if self.eps_tmp is None:
            self.eps_tmp = np.zeros(len(kpt_u[0].eps_n),
                                    dtype=complex)

        # self.Eband = sum_i <psi_i|H|psi_j>
        for kpt in kpt_u:
            self.td_hamiltonian.apply(kpt, kpt.psit_nG, self.hpsit,
                                      calculate_P_ani=False)
            self.mblas.multi_zdotc(self.eps_tmp, kpt.psit_nG,
                                   self.hpsit, len(kpt_u[0].psit_nG))
            self.eps_tmp *= self.wfs.gd.dv
            kpt.eps_n[:] = self.eps_tmp.real

        self.occupations.calculate_band_energy(self.wfs)

        H = self.td_hamiltonian.hamiltonian

        # Nonlocal
        self.Enlxc = 0.0#xcfunc.get_non_local_energy()
        self.Enlkin = H.xc.get_kinetic_energy_correction()

        # PAW
        self.Ekin = H.Ekin0 + self.occupations.e_band + self.Enlkin
        self.Epot = H.Epot
        self.Eext = H.Eext
        self.Ebar = H.Ebar
        self.Exc = H.Exc + self.Enlxc
        self.Etot = self.Ekin + self.Epot + self.Ebar + self.Exc


    def get_td_energy(self):
        """Calculate the time-dependent total energy"""

        if not self.calculate_energy:
            self.Etot = 0.0
            return 0.0

        self.td_overlap.update(self.wfs)
        self.td_density.update()
        self.td_hamiltonian.update(self.td_density.get_density(),
                                   self.time)
        self.update_eigenvalues()

        return self.Etot


    def set_absorbing_boundary(self, absorbing_boundary):
        self.td_hamiltonian.set_absorbing_boundary(absorbing_boundary)


    # exp(ip.r) psi
    def absorption_kick(self, kick_strength):
        """Delta absorption kick for photoabsorption spectrum.

        Parameters:

        kick_strength: [float, float, float]
            Strength of the kick, e.g., [0.0, 0.0, 1e-3]
        
        """
        if self.rank == 0:
            self.text('Delta kick = ', kick_strength)

        self.kick_strength = np.array(kick_strength)

        abs_kick_hamiltonian = AbsorptionKickHamiltonian(self.wfs, self.atoms,
                                   np.array(kick_strength, float))
        abs_kick = AbsorptionKick(self.wfs, abs_kick_hamiltonian,
                                  self.td_overlap, self.solver,
                                  self.preconditioner, self.wfs.gd, self.timer)
        abs_kick.kick()

        # Kick the classical part, if it is present
        if self.hamiltonian.poisson.get_description() == 'FDTD+TDDFT':
            self.hamiltonian.poisson.set_kick(kick = self.kick_strength)



    def __del__(self):
        """Destructor"""
        GPAW.__del__(self)


from gpaw.mpi import world

# Function for calculating photoabsorption spectrum
def photoabsorption_spectrum(dipole_moment_file, spectrum_file,
                             folding='Gauss', width=0.2123,
                             e_min=0.0, e_max=30.0, delta_e=0.05):
    """Calculates photoabsorption spectrum from the time-dependent
    dipole moment.
    
    Parameters:

    dipole_moment_file: string
        Name of the time-dependent dipole moment file from which
        the specturm is calculated
    spectrum_file: string
        Name of the spectrum file
    folding: 'Gauss' or 'Lorentz'
        Whether to use Gaussian or Lorentzian folding
    width: float
        Width of the Gaussian (sigma) or Lorentzian (Gamma)
        Gaussian =     1/(sigma sqrt(2pi)) exp(-(1/2)(omega/sigma)^2)
        Lonrentzian =  (1/pi) (1/2) Gamma / [omega^2 + ((1/2) Gamma)^2]
    e_min: float
        Minimum energy shown in the spectrum (eV)
    e_max: float
        Maxiumum energy shown in the spectrum (eV)
    delta_e: float
        Energy resolution (eV)
    

    """


#    kick_strength: [float, float, float]
#        Strength of the kick, e.g., [0.0, 0.0, 1e-3]
#    fwhm: float
#        Full width at half maximum for peaks in eV
#    delta_omega: float
#        Energy resolution in eV
#    max_energy: float
#        Maximum excitation energy in eV

    if folding != 'Gauss':
        raise RuntimeError('Error in photoabsorption_spectrum: '
                           'Only Gaussian folding is currently supported.')
    
    
    if world.rank == 0:
        print('Calculating photoabsorption spectrum from file "%s"' \
              % dipole_moment_file)

        f_file = file(spectrum_file, 'w')
        dm_file = file(dipole_moment_file, 'r')
        lines = dm_file.readlines()
        dm_file.close()
        # Read kick strength
        columns = lines[0].split('[')
        columns = columns[1].split(']')
        columns = columns[0].split(',')
        kick_strength = np.array([eval(columns[0]),eval(columns[1]),eval(columns[2])], dtype=float)
        strength = np.array(kick_strength, dtype=float)
        # Remove first two lines
        lines.pop(0)
        lines.pop(0)
        print('Using kick strength = ', strength)
        # Continue with dipole moment data
        n = len(lines)
        dm = np.zeros((n,3),dtype=float)
        time = np.zeros((n,),dtype=float)
        for i in range(n):
            data = lines[i].split()
            time[i] = float(data[0])
            dm[i,0] = float(data[2])
            dm[i,1] = float(data[3])
            dm[i,2] = float(data[4])

        t = time - time[0]
        dt = time[1] - time[0]
        dm[:] = dm - dm[0]
        nw = int(e_max / delta_e)
        dw = delta_e * eV_to_aufrequency
        # f(w) = Nw exp(-w^2/2sigma^2)
        #sigma = fwhm / Hartree / (2.* np.sqrt(2.* np.log(2.0)))
        # f(t) = Nt exp(-t^2*sigma^2/2)
        sigma = width * eV_to_aufrequency
        fwhm = sigma * (2.* np.sqrt(2.* np.log(2.0)))
        kick_magnitude = np.sum(strength**2)

        # write comment line
        f_file.write('# Photoabsorption spectrum from real-time propagation\n')
        f_file.write('# GPAW version: ' + str(version) + '\n')
        f_file.write('# Total time = %lf fs, Time step = %lf as\n' \
            % (n * dt * autime_to_attosec/1000.0, \
               dt * autime_to_attosec))
        f_file.write('# Kick = [%lf,%lf,%lf]\n' % (kick_strength[0], \
                                                   kick_strength[1], \
                                                   kick_strength[2]))
        f_file.write('# %sian folding, Width = %lf eV = %lf Hartree <=> ' \
                     'FWHM = %lf eV\n' % (folding, sigma*aufrequency_to_eV, \
                                          sigma, fwhm*aufrequency_to_eV))

        f_file.write('#  om (eV) %14s%20s%20s\n' % ('Sx', 'Sy', 'Sz'))
        # alpha = 2/(2*pi) / eps int dt sin(omega t) exp(-t^2*sigma^2/2)
        #                                * ( dm(t) - dm(0) )
        alpha = 0
        for i in range(nw):
            w = i * dw
            # x
            alphax = np.sum(np.sin(t*w) * np.exp(-t**2*sigma**2/2.0) * dm[:,0])
            alphax *= 2 * dt / (2*np.pi) / kick_magnitude * strength[0]
            # y
            alphay = np.sum(np.sin(t*w) * np.exp(-t**2*sigma**2/2.0) * dm[:,1])
            alphay *= 2 * dt / (2*np.pi) / kick_magnitude * strength[1]
            # z
            alphaz = np.sum(np.sin(t*w) * np.exp(-t**2*sigma**2/2.0) * dm[:,2])
            alphaz *= 2 * dt / (2*np.pi) / kick_magnitude * strength[2]

            # f = 2 * omega * alpha
            line = '%10.6lf %20.10le %20.10le %20.10le\n' \
                % (w*aufrequency_to_eV,
                   2*w*alphax / Hartree,
                   2*w*alphay / Hartree,
                   2*w*alphaz / Hartree)
            f_file.write(line)

            if (i % 100) == 0:
                print('.', end=' ')
                sys.stdout.flush()
                
        print("Sinc contamination", np.exp(-t[-1]**2*sigma**2/2.0))

        print('')
        f_file.close()
        
        print('Calculated photoabsorption spectrum saved to file "%s"' \
              % spectrum_file)

            
    # Make static method
    # photoabsorption_spectrum=staticmethod(photoabsorption_spectrum)
        
