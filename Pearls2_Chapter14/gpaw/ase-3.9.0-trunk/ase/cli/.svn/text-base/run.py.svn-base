from __future__ import division, print_function

import sys
import os.path
import optparse
import os
import tempfile
import time
import traceback
    
import numpy as np

from ase.io import read
from ase.parallel import world
from ase.utils import devnull
from ase.constraints import FixAtoms, UnitCellFilter
from ase.optimize import LBFGS
from ase.io.trajectory import PickleTrajectory
from ase.utils.eos import EquationOfState
from ase.calculators.calculator import get_calculator, names as calcnames
import ase.db as db


def main():
    runner = Runner()
    runner.parse()
    if runner.errors:
        sys.exit(runner.errors)


class Runner:
    def __init__(self):
        self.db = None
        self.opts = None
        self.errors = 0
        self.names = []
        self.calculator_name = None

        if world.rank == 0:
            self.logfile = sys.stdout
        else:
            self.logfile = devnull

    def parse(self, args=None):
        parser = self.make_parser()
        self.add_options(parser)
        self.opts, names = parser.parse_args(args)

        if args is None and self.opts.interactive_python_session:
            file = tempfile.NamedTemporaryFile()
            file.write('import os\n')
            file.write('if "PYTHONSTARTUP" in os.environ:\n')
            file.write('    execfile(os.environ["PYTHONSTARTUP"])\n')
            file.write('from ase.cli.run import Runner\n')
            file.write('atoms = Runner().parse(%r)\n' %
                       ([self.calculator_name] + sys.argv[1:]))
            file.flush()
            os.system('python -i %s' % file.name)
            return

        if self.calculator_name is None:
            if names:
                self.calculator_name = names.pop(0)
            else:
                parser.error('Missing calculator name')
                
        atoms = self.run(names)
        return atoms
        
    def make_parser(self):
        parser = optparse.OptionParser(
            usage='ase-run calculator [options] [system, ...]',
            description="Run calculation with one of ASE's calculators: " +
            ', '.join(calcnames) + '.')
        return parser
        
    def add_options(self, parser):
        add = parser.add_option
        add('-t', '--tag',
            help='String tag added to filenames.')
        add('-p', '--parameters', default='',
            metavar='key=value,...',
            help='Comma-separated key=value pairs of ' +
            'calculator specific parameters.')
        add('-d', '--database',
            help='Use a filename with a ".db" extension for a sqlite3 ' +
            'database or a ".json" extension for a simple json database.  ' +
            'Default is no database')
        add('-S', '--skip', action='store_true',
            help='Skip calculations already done.')
        add('--properties', default='efsdMm',
            help='Default value is "efsdMm" meaning calculate energy, ' +
            'forces, stress, dipole moment, total magnetic moment and ' +
            'atomic magnetic moments.')
        add('-f', '--maximum-force', type=float,
            help='Relax internal coordinates.')
        add('--constrain-tags',
            metavar='T1,T2,...',
            help='Constrain atoms with tags T1, T2, ...')
        add('-s', '--maximum-stress', type=float,
            help='Relax unit-cell and internal coordinates.')
        add('-E', '--equation-of-state', help='Equation of state ...')
        add('--eos-type', default='sjeos', help='Selects the type of eos.')
        add('-i', '--interactive-python-session', action='store_true')
        add('-c', '--collection')
        add('--modify', metavar='...',
            help='Modify atoms with Python statement.  ' +
            'Example: --modify="atoms.positions[-1,2]+=0.1".')
        add('--after', help='Perform operation after calculation.  ' +
            'Example: --after="atoms.calc.write(...)"')

    def log(self, *args, **kwargs):
        print(file=self.logfile, *args, **kwargs)

    def run(self, names):
        opts = self.opts

        if self.db is None:
            # Create database connection:
            self.db = db.connect(opts.database, use_lock_file=True)

        self.expand(names)

        if not names:
            names.insert(0, '-')

        atoms = None
        for name in names:
            if atoms is not None:
                del atoms.calc  # release resources from last calculation
            atoms = self.build(name)
            if opts.modify:
                exec opts.modify in {'atoms': atoms, 'np': np}

            if name == '-':
                name = atoms.info['key_value_pairs']['name']

            skip = False
            id = None
            
            if opts.skip:
                id = self.db.reserve(name=name)
                if id is None:
                    skip = True
        
            if not skip:
                self.set_calculator(atoms, name)

                tstart = time.time()
                try:
                    self.log('Running:', name)
                    data = self.calculate(atoms, name)
                except KeyboardInterrupt:
                    raise
                except Exception:
                    self.log(name, 'FAILED')
                    traceback.print_exc(file=self.logfile)
                    tstop = time.time()
                    data = {'time': tstop - tstart}
                    self.errors += 1
                else:
                    tstop = time.time()
                    data['time'] = tstop - tstart
                    self.db.write(atoms, name=name, data=data)
                
                if id:
                    del self.db[id]

        return atoms
    
    def calculate(self, atoms, name):
        opts = self.opts

        data = {}
        if opts.maximum_force or opts.maximum_stress:
            data = self.optimize(atoms, name)
        if opts.equation_of_state:
            data.update(self.eos(atoms, name))
        data.update(self.calculate_once(atoms, name))

        if opts.after:
            exec opts.after in {'atoms': atoms, 'data': data}

        return data

    def expand(self, names):
        if not self.names and self.opts.collection:
            con = db.connect(self.opts.collection)
            self.names = [dct.id for dct in con.select()]
        if not names:
            names[:] = self.names
            return
        if not self.names:
            return
        i = 0
        while i < len(names):
            name = names[i]
            if name.count('-') == 1:
                s1, s2 = name.split('-')
                if s1 in self.names and s2 in self.names:
                    j1 = self.names.index(s1)
                    j2 = self.names.index(s2)
                    names[i:i + 1] = self.names[j1:j2 + 1]
                    i += j2 - j1
            i += 1

    def build(self, name):
        if name == '-':
            con = db.connect(sys.stdin, 'json')
            return con.get_atoms(add_additional_information=True)
        elif self.opts.collection:
            con = db.connect(self.opts.collection)
            return con.get_atoms(name)
        else:
            return read(name)

    def set_calculator(self, atoms, name):
        cls = get_calculator(self.calculator_name)
        parameters = str2dict(self.opts.parameters)
        if getattr(cls, 'nolabel', False):
            atoms.calc = cls(**parameters)
        else:
            atoms.calc = cls(label=self.get_filename(name), **parameters)

    def calculate_once(self, atoms, name):
        opts = self.opts

        for p in opts.properties or 'efsdMm':
            property, method = {'e': ('energy', 'get_potential_energy'),
                                'f': ('forces', 'get_forces'),
                                's': ('stress', 'get_stress'),
                                'd': ('dipole', 'get_dipole_moment'),
                                'M': ('magmom', 'get_magnetic_moment'),
                                'm': ('magmoms', 'get_magnetic_moments')}[p]
            try:
                getattr(atoms, method)()
            except NotImplementedError:
                pass

        data = {}
        
        return data

    def optimize(self, atoms, name):
        opts = self.opts
        if opts.constrain_tags:
            tags = [int(t) for t in opts.constrain_tags.split(',')]
            mask = [t in tags for t in atoms.get_tags()]
            atoms.constraints = FixAtoms(mask=mask)
        
        trajectory = PickleTrajectory(self.get_filename(name, 'traj'), 'w',
                                      atoms)
        if opts.maximum_stress:
            optimizer = LBFGS(UnitCellFilter(atoms), logfile=self.logfile)
            fmax = opts.maximum_stress
        else:
            optimizer = LBFGS(atoms, logfile=self.logfile)
            fmax = opts.maximum_force

        optimizer.attach(trajectory)
        optimizer.run(fmax=fmax)

        data = {}
        if hasattr(optimizer, 'force_calls'):
            data['force_calls'] = optimizer.force_calls

        return data

    def eos(self, atoms, name):
        opts = self.opts
        
        traj = PickleTrajectory(self.get_filename(name, 'traj'), 'w', atoms)
        eps = 0.01
        strains = np.linspace(1 - eps, 1 + eps, 5)
        v1 = atoms.get_volume()
        volumes = strains**3 * v1
        energies = []
        cell1 = atoms.cell
        for s in strains:
            atoms.set_cell(cell1 * s, scale_atoms=True)
            energies.append(atoms.get_potential_energy())
            traj.write(atoms)
        traj.close()
        eos = EquationOfState(volumes, energies, opts.eos_type)
        v0, e0, B = eos.fit()
        atoms.set_cell(cell1 * (v0 / v1)**(1 / 3), scale_atoms=True)
        data = {'volumes': volumes,
                'energies': energies,
                'fitted_energy': e0,
                'fitted_volume': v0,
                'bulk_modulus': B,
                'eos_type': opts.eos_type}
        return data

    def get_filename(self, name=None, ext=None):
        if name is None:
            if self.opts.tag is None:
                filename = 'ase'
            else:
                filename = self.opts.tag
        else:
            if '.' in name:
                name = name.rsplit('.', 1)[0]
            if self.opts.tag is None:
                filename = name
            else:
                filename = name + '-' + self.opts.tag

        if ext:
            filename += '.' + ext

        return filename


def str2dict(s, namespace={}, sep='='):
    """Convert comma-separated key=value string to dictionary.

    Examples:

    >>> str2dict('xc=PBE,nbands=200,parallel={band:4}')
    {'xc': 'PBE', 'nbands': 200, 'parallel': {'band': 4}}
    >>> str2dict('a=1.2,b=True,c=ab,d=1,2,3,e={f:42,g:cd}')
    {'a': 1.2, 'c': 'ab', 'b': True, 'e': {'g': 'cd', 'f': 42}, 'd': (1, 2, 3)}
    """
    
    def myeval(value):
        try:
            value = eval(value, namespace)
        except (NameError, SyntaxError):
            pass
        return value

    dct = {}
    s = (s + ',').split(sep)
    for i in range(len(s) - 1):
        key = s[i]
        m = s[i + 1].rfind(',')
        value = s[i + 1][:m]
        if value[0] == '{':
            assert value[-1] == '}'
            value = str2dict(value[1:-1], namespace, ':')
        elif value[0] == '(':
            assert value[-1] == ')'
            value = [myeval(t) for t in value[1:-1].split(',')]
        else:
            value = myeval(value)
        dct[key] = value
        s[i + 1] = s[i + 1][m + 1:]
    return dct
