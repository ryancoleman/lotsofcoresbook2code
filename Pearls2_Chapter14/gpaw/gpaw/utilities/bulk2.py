"""This module contains helper classes for running simple calculations.

XXX this file should be renamed to something else!!!

The main user of the Runner classes defined in this module is the
``gpaw`` command-line tool.  """

import os
import sys
from math import sqrt

import numpy as np
from ase.atoms import Atoms, string2symbols
from ase.utils.eos import EquationOfState
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io.trajectory import PickleTrajectory
from ase.optimize.lbfgs import LBFGS
from ase.constraints import FixAtoms
from ase.io import read
import ase.units as units
from ase.data import covalent_radii
from ase.utils import devnull

from gpaw.aseinterface import GPAW
from gpaw.poisson import PoissonSolver
from gpaw.mpi import world
from gpaw.utilities import h2gpts
from gpaw.occupations import FermiDirac

# Magnetic moments of isolated atoms:
magmom = {'C': 2, 'N': 3, 'Pt': 2, 'F': 1, 'Mg': 0, 'Na': 1, 'Cl': 1, 'Al': 1,
          'O': 2, 'Li': 1, 'P': 3, 'Si': 2, 'Cu': 1, 'Fe': 4}

class Runner:
    """Base class for running the calculations.

    Subclasses must implement a set_calculator() method."""
    
    def __init__(self, name, atoms, strains=None, tag='', clean=True,
                 fmax=None, out='-'):
        """Construct runner object.

        Results will be written to trajectory files or read from those
        files if the already exist.
        
        name: str
            Name of calculation.
        atoms: Atoms object
            configuration to work on.
        strains: list of floats
            The list of strains to apply to the unit cell or bond length.
            Defaults to [1.0].
        tag: str
            Tag used for filenames like <name>-<tag>.traj.
        clean: bool
            Do *not* read results from files.
        fmax: float
            Optimize geometry until all forces are below fmax.
        """
        
        if strains is None:
            strains = [1.0]

        if tag:
            self.tag = '-' + tag
        else:
            self.tag = ''

        if world.rank == 0:
            if out is None:
                out = devnull
            elif isinstance(out, str):
                if out == '-':
                    out = sys.stdout
                else:
                    out = open(out, 'w')
        else:
            out = devnull
            
        self.name = name
        self.atoms = atoms
        self.strains = np.array(strains)
        self.clean = clean
        self.out = out
        self.fmax = fmax
        
        self.volumes = None
        self.bondlengths = None
        self.energies = None

        self.dimer = (not atoms.pbc.any() and len(self.atoms) == 2)

    def log(self, *args, **kwargs):
        self.out.write(kwargs.get('sep', ' ').join([str(arg)
                                                    for arg in args]) +
                       kwargs.get('end', '\n'))

    def run(self):
        """Start calculation or read results from file."""
        filename = '%s%s.traj' % (self.name, self.tag)
        if self.clean or not os.path.isfile(filename):
            world.barrier()
            if world.rank == 0:
                open(filename, 'w')
            self.calculate(filename)
        else:
            try:
                self.log('Reading', filename, end=' ')
                configs = read(filename, ':')
            except IOError:
                self.log('FAILED')
            else:
                self.log()
                if self.fmax is not None:
                    configs = configs[-1:]
                # Extract bondlengths/volumes and energies:
                if self.dimer:
                    self.bondlengths = [a.get_distance(0, 1) for a in configs]
                else:
                    self.volumes = [a.get_volume() for a in configs]
                self.energies = [a.get_potential_energy() for a in configs]

    def calculate(self, filename):
        """Run calculation and write results to file."""
        self.log('Calculating', self.name, '...')
        config = self.atoms.copy()
        self.set_calculator(config, filename)
        traj = PickleTrajectory(filename, 'w', backup=False)
        cell = config.get_cell()
        self.energies = []
        if self.fmax is not None:
            if (not config.constraints and
                not config.positions.any(axis=1).all()):
                # one atom is at (0,0,0) - fix it:
                mask = np.logical_not(config.positions.any(axis=1))
                config.constraints = FixAtoms(mask=mask)
            dyn = LBFGS(config)
            dyn.attach(traj, 1, config)
            dyn.run(fmax=self.fmax)
            e = config.get_potential_energy()
            self.post_process(config)
            self.energies.append(e)
        elif not config.pbc.any() and len(config) == 2:
            # This is a dimer.
            self.bondlengths = []
            d0 = config.get_distance(0, 1)
            for strain in self.strains:
                d = d0 * strain
                config.set_distance(0, 1, d)
                self.bondlengths.append(d)
                e = config.get_potential_energy()
                self.energies.append(e)
                self.post_process(config)
                traj.write(config)
        else:
            self.volumes = []
            for strain in self.strains:
                config.set_cell(strain * cell, scale_atoms=True)
                self.volumes.append(config.get_volume())
                e = config.get_potential_energy()
                self.energies.append(e)
                self.post_process(config)
                traj.write(config)
        return config
    
    def post_process(self, config):
        pass

    def summary(self, plot=False, a0=None):
        if self.energies is None:
            return
        natoms = len(self.atoms)
        if len(self.energies) == 1:
            self.log('Total energy: %.3f eV (%d atom%s)' %
                     (self.energies[0], natoms, ' s'[1:natoms]))
        elif self.dimer:
            return self.dimer_summary(plot)
        else:
            return self.bulk_summary(plot, a0)

    def dimer_summary(self, plot):
        d = np.array(self.bondlengths)
        fit0 = np.poly1d(np.polyfit(1 / d, self.energies, 3))
        fit1 = np.polyder(fit0, 1)
        fit2 = np.polyder(fit1, 1)

        d0 = None
        for t in np.roots(fit1):
            if t > 0 and fit2(t) > 0:
                d0 = 1 / t
                break

        if d0 is None:
            raise ValueError('No minimum!')
        
        e0 = fit0(t)
        k = fit2(t) * t**4
        m1, m2 = self.atoms.get_masses()
        m = m1 * m2 / (m1 + m2)
        hnu = units._hbar * 1e10 * sqrt(k / units._e / units._amu / m)

        self.log('Fit using %d points:' % len(self.energies))
        self.log('Bond length: %.3f Ang^3' % d0)
        self.log('Frequency: %.1f meV (%.1f cm^-1)' %
                 (1000 * hnu,
                  hnu * 0.01 * units._e / units._c / units._hplanck))
        self.log('Total energy: %.3f eV' % e0)

        if plot:
            import pylab as plt
            plt.plot(d, self.energies, 'o')
            x = np.linspace(d[0], d[-1], 50)
            plt.plot(x, fit0(x**-1), '-r')
            plt.show()

        dimer = self.atoms.copy()
        dimer.set_distance(0, 1, d0)
        self.write_optimized(dimer, e0)
            
        return e0, d0, hnu

    def bulk_summary(self, plot, a0):
        natoms = len(self.atoms)
        eos = EquationOfState(self.volumes, self.energies)
        v, e, B = eos.fit()
        x = (v / self.atoms.get_volume())**(1.0 / 3)

        self.log('Fit using %d points:' % len(self.energies))
        self.log('Volume per atom: %.3f Ang^3' % (v / natoms))
        if a0:
            a = a0 * x
            self.log('Lattice constant: %.3f Ang' % a)
        else:
            a = None
        self.log('Bulk modulus: %.1f GPa' % (B * 1e24 / units.kJ))
        self.log('Total energy: %.3f eV (%d atom%s)' %
                 (e, natoms, ' s'[1:natoms]))

        if plot:
            import pylab as plt
            plt.plot(self.volumes, self.energies, 'o')
            x = np.linspace(self.volumes[0], self.volumes[-1], 50)
            plt.plot(x, eos.fit0(x**-(1.0 / 3)), '-r')
            plt.show()
            
        bulk = self.atoms.copy()
        bulk.set_cell(x * bulk.cell, scale_atoms=True)
        self.write_optimized(bulk, e)

        return e, v, B, a

    def write_optimized(self, atoms, energy):
        magmoms = None  # XXX we could do better ...
        forces = np.zeros((len(atoms), 3))
        calc = SinglePointCalculator(energy, forces, np.zeros((3, 3)),
                                     magmoms, atoms)
        atoms.set_calculator(calc)
        filename = '%s%s-optimized.traj' % (self.name, self.tag)
        traj = PickleTrajectory(filename, 'w', backup=False)
        traj.write(atoms)
        traj.close()


class EMTRunner(Runner):
    """EMT implementation"""
    def set_calculator(self, config, filename):
        config.set_calculator(EMT())


class GPAWRunner(Runner):
    """GPAW implementation"""
    def set_parameters(self, vacuum=3.0, write_gpw_file=None, **kwargs):
        self.vacuum = vacuum
        self.write_gpw_file = write_gpw_file
        self.gpwfilename = None
        self.input_parameters = kwargs
        
    def set_calculator(self, config, filename):
        kwargs = {}
        kwargs.update(self.input_parameters)

        if self.write_gpw_file is not None:
            self.gpwfilename = filename[:-4] + 'gpw'

        if 'txt' not in kwargs:
            kwargs['txt'] = filename[:-4] + 'txt'
        
        if not config.pbc.any():
            # Isolated atom or molecule:
            config.center(vacuum=self.vacuum)
            if (len(config) == 1 and
                config.get_initial_magnetic_moments().any()):
                kwargs['hund'] = True

        # Use fixed number of gpts:
        if 'gpts' not in kwargs:
            h = kwargs.get('h', 0.2)
            gpts = h2gpts(h, config.cell)
            kwargs['h'] = None
            kwargs['gpts'] = gpts
        
        self.calc = GPAW(**kwargs)
        config.set_calculator(self.calc)

    def get_calculator(self):
        """Returns the calculator object - available when finished only."""
        return self.calc

    def post_process(self, config):
        if self.write_gpw_file is not None:
            self.calc.write(self.gpwfilename, mode=self.write_gpw_file)
        
    def check_occupation_numbers(self, config):
        """Check that occupation numbers are integers."""
        if config.pbc.any():
            return
        calc = config.get_calculator()
        nspins = calc.get_number_of_spins()
        for s in range(nspins):
            f = calc.get_occupation_numbers(spin=s)
            if abs(f % (2 // nspins)).max() > 0.0001:
                raise RuntimeError('Fractional occupation numbers?!')
