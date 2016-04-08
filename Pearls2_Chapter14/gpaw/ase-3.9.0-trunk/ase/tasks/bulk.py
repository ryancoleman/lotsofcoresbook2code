import optparse

import numpy as np

from ase.lattice import bulk
from ase.constraints import StrainFilter
import ase.optimize
from ase.tasks.task import OptimizeTask
from ase.data import chemical_symbols, reference_states
from ase.utils.eos import EquationOfState
from ase.io.trajectory import PickleTrajectory


class BulkTask(OptimizeTask):
    taskname = 'bulk'

    def __init__(self, crystal_structure=None, lattice_constant=None,
                 c_over_a=None, cubic=False, orthorhombic=False, fit=None,
                 eos=None, sfmax=None, soptimizer='BFGS', ssteps=100000000,
                 **kwargs):
        """Bulk task."""

        self.crystal_structure = crystal_structure
        self.lattice_constant = lattice_constant
        self.c_over_a = c_over_a
        self.cubic = cubic
        self.orthorhombic = orthorhombic
        self.eos = eos
        self.fit = fit
        self.sfmax = sfmax
        self.soptimizer = soptimizer
        self.ssteps = ssteps

        self.repeat = None

        OptimizeTask.__init__(self, **kwargs)

        self.summary_keys = ['energy', 'fitted energy', 'volume', 'B']

    def expand(self, names):
        """Expand fcc, bcc, hcp and diamond.

        The name fcc will be expanded to all the elements with the fcc
        stucture and so on."""

        names = OptimizeTask.expand(self, names)

        newnames = []
        for name in names:
            if name in ['fcc', 'bcc', 'hcp', 'diamond']:
                for Z in range(1, 95):
                    x = reference_states[Z]
                    if x is not None and x['symmetry'] == name:
                        newnames.append(chemical_symbols[Z])
            else:
                newnames.append(name)

        return newnames

    def build_system(self, name):
        atoms = bulk(name, crystalstructure=self.crystal_structure,
                     a=self.lattice_constant, covera=self.c_over_a,
                     orthorhombic=self.orthorhombic, cubic=self.cubic)

        M = {'Fe': 2.3, 'Co': 1.2, 'Ni': 0.6}.get(name)
        if M is not None:
            atoms.set_initial_magnetic_moments([M] * len(atoms))

        if self.repeat is not None:
            r = self.repeat.split(',')
            if len(r) == 1:
                r = 3 * r
            atoms = atoms.repeat([int(c) for c in r])

        return atoms

    def fit_volume(self, name, atoms, data=None):
        N, x = self.fit
        cell0 = atoms.get_cell()
        v = atoms.get_volume()
        if x > 0:
            strains = np.linspace(1 - x, 1 + x, N)
        else:
            strains = np.linspace(1 + x / v, 1 - x / v,  N)**(1./3)
        energies = []
        traj = PickleTrajectory(self.get_filename(name, 'fit.traj'), 'w')
        for s in strains:
            atoms.set_cell(cell0 * s, scale_atoms=True)
            energies.append(atoms.get_potential_energy())
            traj.write(atoms)

        traj.close()

        if data is not None:
            data['strains'] = strains
            data['energies'] = energies
        else:
            assert N % 2 == 1
            data = {'energy': energies[N // 2],
                    'strains': strains,
                    'energies': energies}

        return data

    def soptimize(self, name, atoms, data, trajectory=None):
        # Call it soptimize to avoid conflicts with OptimizeTask.optimize.
        # With so many levels of inheritance one never knows ...
        optstr = "ase.optimize." + self.soptimizer
        optstr += "(atoms, logfile=self.logfile)"
        optimizer = eval(optstr)
        if trajectory is not None:
            if not isinstance(trajectory, str):
                optimizer.attach(trajectory)
        try:
            # handle scipy optimizers who raise Converged when done
            from ase.optimize import Converged
            try:
                optimizer.run(fmax=self.sfmax, steps=self.ssteps)
            except Converged:
                raise
        except ImportError:
            optimizer.run(fmax=self.sfmax, steps=self.ssteps)
        # StrainFilter optimizer steps
        steps = optimizer.get_number_of_steps()
        if data.get('strain optimizer steps', None) is None:
            data['strain optimizer steps'] = steps
        else:
            data['strain optimizer steps'] += steps
        # optimizer force calls
        if hasattr(optimizer, 'force_calls'):
            calls = optimizer.force_calls
        else:
            calls = steps
        if data.get('strain optimizer force calls', None) is None:
            data['strain optimizer force calls'] = calls
        else:
            data['strain optimizer force calls'] += calls

    def converged(self, atoms, sfmax, fmax):
        # The same criteria as in ASE optimizers:
        # see ase.constraints: StrainFilter
        sforces = - atoms.get_stress().ravel()*atoms.get_volume()
        rmssforces = np.sum(sforces**2)**0.5
        maxsforce = rmssforces
        # see ase.optimize.optimize: converged
        rmsforces = np.sum(atoms.get_forces()**2, axis=1)**0.5
        maxforce = max(rmsforces)

        if maxforce < fmax and maxsforce < sfmax:
            return True
        return False

    def calculate(self, name, atoms):
        #????
        if self.sfmax is not None and self.fmax is not None:
            # this performs first relaxation of internal degrees of freedom
            data = OptimizeTask.calculate(self, name, atoms)
            # writing traj from optimizer does not work for StrainFilter!
            traj = PickleTrajectory(self.get_filename(name, 'traj'), 'a', atoms)
            sf = StrainFilter(atoms)
            while not self.converged(atoms, sfmax=self.sfmax, fmax=self.fmax):
                # take a step on the cell
                self.soptimize(name, sf, data, trajectory=traj)
                # relax internal degrees of freedom
                OptimizeTask.optimize(self, name, atoms, data, trajectory=traj)
            data['relaxed energy'] = atoms.get_potential_energy()
            data['relaxed volume'] = atoms.get_volume()
        elif self.sfmax is not None:
            # this performs single-point energy calculation
            data = OptimizeTask.calculate(self, name, atoms)
            sf = StrainFilter(atoms)
            # writing traj from optimizer does not work for StrainFilter!
            traj = PickleTrajectory(self.get_filename(name, 'traj'), 'w', atoms)
            self.soptimize(name, sf, data, trajectory=traj)
            data['relaxed energy'] = atoms.get_potential_energy()
            data['relaxed volume'] = atoms.get_volume()
        elif self.fmax is not None:
            data = OptimizeTask.calculate(self, name, atoms)
        else:
            # no optimization
            if self.fit is None:
                # only calculate single-point energy if no fit follows
                data = OptimizeTask.calculate(self, name, atoms)
        if self.fit is not None:
            if self.sfmax is not None or self.fmax is not None:
                # fit after optimization
                self.fit_volume(name, atoms, data)
            else:
                # fit is the only task performed
                data = self.fit_volume(name, atoms)
        return data

    def analyse(self):
        for name, data in self.data.items():
            if 'strains' in data:
                atoms = self.create_system(name)
                # use relaxed volume if present
                if 'relaxed volume' in data:
                    volume = data['relaxed volume']
                else:
                    volume = atoms.get_volume()
                volumes = data['strains']**3 * volume
                energies = data['energies']
                # allow selection of eos type independent of data
                if self.eos is not None:
                    eos = EquationOfState(volumes, energies, self.eos)
                else:
                    eos = EquationOfState(volumes, energies)
                try:
                    v, e, B = eos.fit()
                except (RuntimeError, ValueError):
                    pass
                else:
                    data['fitted energy'] = e
                    data['volume'] = v
                    data['B'] = B

                    if abs(v) < min(volumes) or abs(v) > max(volumes):
                        raise ValueError(name + ': fit outside of range! ' + \
                                         str(abs(v)) + ' not in ' + \
                                         str(volumes))

    def add_options(self, parser):
        OptimizeTask.add_options(self, parser)

        bulk = optparse.OptionGroup(parser, 'Bulk')
        bulk.add_option('-F', '--fit', metavar='N,x',
                        help='Find optimal volume and bulk modulus ' +
                        'using odd N points and variations of the lattice ' +
                        'constant a from -x % to +x %, i.e. in the interval '
                        '<a - a * x * 100, ..., a, ..., a + a * x * 100>. ' +
                        'This method gives non-equidistant sampling of volume. ' +
                        'With x negative (in Angstrom**3) the sampling of ' +
                        'the cell volume (v) in the interval ' +
                        '<(1 + x /v), ..., 1, ..., (1 - x /v)> is used. ' +
                        'This method gives equidistant sampling of volume.')
        bulk.add_option('--eos', type='str',
                        metavar='eos',
                        help='Selects the type of eos.')
        bulk.add_option('--srelax', metavar='SFMAX[,SOPTIMIZER]',
                        help='Relax cell by minimizing stress using StranFilter '
                        'with SOPTIMIZER algorithm. The SOPTIMIZER keyword is '
                        'optional, and if omitted BFGS is used by default.')
        bulk.add_option('--srelaxsteps', type='int',
                        metavar='ssteps',
                        help='Limit the number of SF optimizer steps.')
        bulk.add_option('-x', '--crystal-structure',
                        help='Crystal structure.',
                        choices=['sc', 'fcc', 'bcc', 'hcp', 'diamond',
                                 'zincblende', 'rocksalt', 'cesiumchloride',
                                 'fluorite'])
        bulk.add_option('-a', '--lattice-constant', type='float',
                        help='Lattice constant in Angstrom.')
        bulk.add_option('--c-over-a', type='float',
                        help='c/a ratio.')
        bulk.add_option('-O', '--orthorhombic', action='store_true',
                        help='Use orthorhombic unit cell.')
        bulk.add_option('-C', '--cubic', action='store_true',
                        help='Use cubic unit cell.')
        bulk.add_option('-r', '--repeat',
                        help='Repeat unit cell.  Use "-r 2" or "-r 2,3,1".')
        parser.add_option_group(bulk)

    def parse(self, opts, args):
        OptimizeTask.parse(self, opts, args)

        if opts.srelax:
            if len(opts.srelax.split(',')) > 1:
                self.sfmax, self.soptimizer = opts.srelax.split(',')
            else:
                self.sfmax = opts.srelax
                self.soptimizer = 'BFGS'
            self.sfmax = float(self.sfmax)

        if opts.srelaxsteps is not None:
            self.ssteps = int(opts.srelaxsteps)
        else:
            # yes, the default number of ASE optimizer steps
            # ase/optimize/optimize.py
            self.ssteps = 100000000

        if opts.fit:
            points, strain = opts.fit.split(',')
            if float(strain) > 0:
                self.fit = (int(points), float(strain) * 0.01)
            else:
                self.fit = (int(points), float(strain))

        self.eos = opts.eos

        self.crystal_structure = opts.crystal_structure
        self.lattice_constant = opts.lattice_constant
        self.c_over_a = opts.c_over_a
        self.orthorhombic = opts.orthorhombic
        self.cubic = opts.cubic
        self.repeat = opts.repeat
