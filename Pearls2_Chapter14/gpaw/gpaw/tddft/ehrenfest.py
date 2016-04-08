from ase.units import Bohr, AUT, _me, _amu
from gpaw import *
from gpaw.tddft import *
from gpaw.tddft.units import attosec_to_autime
from gpaw.mpi import world

###############################################################################
# EHRENFEST DYNAMICS WITHIN THE PAW METHOD
# WORKS WITH PAW AS LONG THERE ISN'T TOO
# TOO MUCH OVERLAP BETWEEN THE SPHERES
# SUPPORTS ALSO HGH PSEUDOPOTENTIALS
###############################################################################

# m a(t+dt)   = F[psi(t),x(t)] 
# x(t+dt/2)   = x(t) + v(t) dt/2 + .5 a(t) (dt/2)^2
# vh(t+dt/2)  = v(t) + .5 a(t) dt/2
# m a(t+dt/2) = F[psi(t),x(t+dt/2)] 
# v(t+dt/2)   = vh(t+dt/2) + .5 a(t+dt/2) dt/2
#
# psi(t+dt)   = U(t,t+dt) psi(t)
#
# m a(t+dt/2) = F[psi(t+dt),x(t+dt/2)] 
# x(t+dt)     = x(t+dt/2) + v(t+dt/2) dt/2 + .5 a(t+dt/2) (dt/2)^2
# vh(t+dt)    = v(t+dt/2) + .5 a(t+dt/2) dt/2
# m a(t+dt)   = F[psi(t+dt),x(t+dt)] 
# v(t+dt)     = vh(t+dt/2) + .5 a(t+dt/2) dt/2

#TODO: move force corrections from forces.py to this module, as well as
# the cg method for calculating the inverse of S from overlap.py

class EhrenfestVelocityVerlet:
    
    def __init__(self, calc, mass_scale = 1.0, setups='paw'):
        """Initializes the Ehrenfest MD calculator.

        Parameters
        ----------

        calc: TDDFT Object

        mass_scale: 1.0
            Scaling coefficient for atomic masses

        setups: {'paw', 'hgh'}
            Type of setups to use for the calculation

        Note
        ------

        Use propagator = 'EFSICN' for when creating the TDDFT object from a PAW ground state
        calculator and propagator = 'EFSICN_HGH' for HGH pseudopotentials

        """
        #print '--- EhrenfestVelocityVerlet is NOT READY FOR PRODUCTION USE ---'
        self.calc = calc
        self.setups = setups
        self.x  = self.calc.atoms.positions.copy() / Bohr
        self.xn = self.x.copy()
        self.v  = self.x.copy()
        amu_to_aumass = _amu / _me
        if self.calc.atoms.get_velocities() is not None:
            self.v = self.calc.atoms.get_velocities() / (Bohr / AUT)
        else:
            self.v[:] = 0.0
            self.calc.atoms.set_velocities(self.v)
        
        self.vt = self.v.copy()
        self.vh = self.v.copy()
        self.time = 0.0
        
        self.M = calc.atoms.get_masses() * amu_to_aumass * mass_scale

        self.a  = self.v.copy()
        self.ah = self.a.copy()
        self.an = self.a.copy()
        self.F  = self.a.copy()

        self.calc.get_td_energy()
        self.calc.forces.reset()
        self.F = self.calc.forces.calculate(self.calc.wfs,
                                            self.calc.td_density.get_density(),
                                            self.calc.td_hamiltonian.hamiltonian)

        for i in range(len(self.F)):
            self.a[i] = self.F[i] / self.M[i]


    def propagate(self, dt):
        """Performs one Ehrenfest MD propagation step

        Parameters
        ----------

        dt: scalar
            Time step (in attoseconds) used for the Ehrenfest MD step

        """
        self.x  = self.calc.atoms.positions.copy() / Bohr
        self.v  = self.calc.atoms.get_velocities() / (Bohr / AUT)

        dt = dt * attosec_to_autime

        # m a(t+dt)   = F[psi(t),x(t)] 
        self.calc.atoms.positions = self.x * Bohr
        self.calc.set_positions(self.calc.atoms)
        self.calc.get_td_energy()
        self.calc.forces.reset()
        self.F = self.calc.forces.calculate(self.calc.wfs,
                                            self.calc.td_density.get_density(),
                                            self.calc.td_hamiltonian.hamiltonian)

        for i in range(len(self.F)):
            self.a[i] = self.F[i] / self.M[i]

        # x(t+dt/2)   = x(t) + v(t) dt/2 + .5 a(t) (dt/2)^2
        # vh(t+dt/2)  = v(t) + .5 a(t) dt/2
        self.xh  = self.x + self.v * dt/2 + .5 * self.a * dt/2*dt/2
        self.vhh = self.v + .5 * self.a * dt/2

        # m a(t+dt/2) = F[psi(t),x(t+dt/2)a] 
        self.calc.atoms.positions = self.xh * Bohr
        self.calc.set_positions(self.calc.atoms)
        self.calc.get_td_energy()
        self.calc.forces.reset()
        self.F = self.calc.forces.calculate(self.calc.wfs,
                                            self.calc.td_density.get_density(),
                                            self.calc.td_hamiltonian.hamiltonian)

        for i in range(len(self.F)):
            self.ah[i] = self.F[i] / self.M[i]

        # v(t+dt/2)   = vh(t+dt/2) + .5 a(t+dt/2) dt/2
        self.vh = self.vhh + .5 * self.ah * dt/2


        # Propagate wf
        # psi(t+dt)   = U(t,t+dt) psi(t)
        if(self.setups == 'paw'):
            niters = self.calc.propagator.propagate(self.time, dt, self.vh)
        else:
            niters = self.calc.propagator.propagate(self.time, dt)
        #print 'Propagation took = ', niters

        # m a(t+dt/2) = F[psi(t+dt),x(t+dt/2)] 
        self.calc.atoms.positions = self.xh * Bohr
        self.calc.set_positions(self.calc.atoms)
        self.calc.get_td_energy()
        self.calc.forces.reset()
        self.F = self.calc.forces.calculate(self.calc.wfs,
                                            self.calc.td_density.get_density(),
                                            self.calc.td_hamiltonian.hamiltonian)

        for i in range(len(self.F)):
            self.ah[i] = self.F[i] / self.M[i]

        # x(t+dt)     = x(t+dt/2) + v(t+dt/2) dt/2 + .5 a(t+dt/2) (dt/2)^2
        # vh(t+dt)    = v(t+dt/2) + .5 a(t+dt/2) dt/2
        self.xn  = self.xh + self.vh * dt/2 + .5 * self.ah * dt/2*dt/2
        self.vhh = self.vh + .5 * self.ah * dt/2

        # m a(t+dt)   = F[psi(t+dt),x(t+dt)] 
        self.calc.atoms.positions = self.xn * Bohr
        self.calc.set_positions(self.calc.atoms)
        self.calc.get_td_energy()
        self.calc.update_eigenvalues()
        self.calc.forces.reset()
        self.F = self.calc.forces.calculate(self.calc.wfs,
                                            self.calc.td_density.get_density(),
                                            self.calc.td_hamiltonian.hamiltonian)

        for i in range(len(self.F)):
            self.an[i] = self.F[i] / self.M[i]

        # v(t+dt)     = vh(t+dt/2) + .5 a(t+dt/2) dt/2
        self.vn = self.vhh + .5 * self.an * dt/2

        
        #print '--- X ---'
        #print self.x
        #print self.xn
        #print self.xn - self.x
        #print '-- V ---'
        #print self.v
        #print self.vn
        #print self.vn - self.v
        #print '--- A ---'
        #print self.a
        #print self.an
        #print self.an - self.a
        #print '---'

        # update
        self.x[:] = self.xn
        self.v[:] = self.vn
        self.a[:] = self.an

        # update atoms
        self.calc.atoms.set_positions(self.x * Bohr)
        self.calc.atoms.set_velocities(self.v * Bohr / AUT)

    def get_energy(self):
        """Updates kinetic, electronic and total energies"""
        self.Ekin = 0.0
        for i in range(len(self.v)):
            self.Ekin += (
                .5 * self.M[i] * self.v[i][0]**2
                +.5 * self.M[i] * self.v[i][1]**2
                +.5 * self.M[i] * self.v[i][2]**2
                )
        self.Epot = self.calc.get_td_energy()
        self.Etot = self.Ekin + self.Epot

        #print 'Ecur = ', [self.Etot, self.Ekin, self.Epot]
        
        return self.Etot
        
    def get_velocities_in_au(self):
        return self.v

    def set_velocities_in_au(self, v):
        self.v[:] = v
        self.calc.atoms.set_velocities(v * Bohr / AUT)
