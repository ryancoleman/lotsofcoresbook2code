import optparse
from math import sqrt, pi

import numpy as np

from ase.dft.kpoints import monkhorst_pack


def str2dict(s, namespace={}, sep='='):
    """Convert comma-separated key=value string to dictionary.

    Examples:

    >>> str2dict('xc=PBE,nbands=200,parallel={band:4}')
    {'xc': 'PBE', 'nbands': 200, 'parallel': {'band': 4}}
    >>> str2dict('a=1.2,b=True,c=ab,d=1,2,3,e={f:42,g:cd}')
    {'a': 1.2, 'c': 'ab', 'b': True, 'e': {'g': 'cd', 'f': 42}, 'd': (1, 2, 3)}
    """

    dct = {}
    s = (s + ',').split(sep)
    for i in range(len(s) - 1):
        key = s[i]
        m = s[i + 1].rfind(',')
        value = s[i + 1][:m]
        if value[0] == '{':
            assert value[-1] == '}'
            value = str2dict(value[1:-1], {}, ':')
        else:
            try:
                value = eval(value, namespace)
            except (NameError, SyntaxError):
                pass
        dct[key] = value
        s[i + 1] = s[i + 1][m + 1:]
    return dct


class CalculatorFactory:
    def __init__(self, Class, name, label='label',
                 kpts=None, kptdensity=3.0,
                 **kwargs):
        """Calculator factory object.

        Used to create calculators with specific parameters."""

        self.Class = Class
        self.name = name
        self.label = label
        self.kpts = kpts
        self.kptdensity = kptdensity
        self.kwargs = kwargs

    def calculate_kpts(self, atoms):
        """Estimate an appropriate number of k-points."""

        if self.kpts is not None:
            # Number of k-points was explicitely set:
            return self.kpts

        # Use kptdensity to make a good estimate:
        recipcell = atoms.get_reciprocal_cell()
        kpts = []
        for i in range(3):
            if atoms.pbc[i]:
                k = 2 * pi * sqrt((recipcell[i]**2).sum()) * self.kptdensity
                kpts.append(max(1, 2 * int(round(k / 2))))
            else:
                kpts.append(1)

        return kpts

    def __call__(self, name, atoms):
        """Create calculator.

        Put name in the filename of all created files."""

        kpts = self.calculate_kpts(atoms)
        if kpts != 'no k-points':
            if self.name == 'aims':  # XXX Aims uses k_grid!
                self.kwargs['k_grid'] = kpts
            else:
                self.kwargs['kpts'] = kpts

        if self.label is not None:
            self.kwargs[self.label] = name
        return self.Class(**self.kwargs)
        
        if self.label is None:
            return self.Class(**self.kwargs)
        else:
            return self.Class(name, **self.kwargs)

    def add_options(self, parser):
        calc = optparse.OptionGroup(parser, 'Calculator')
        calc.add_option('-k', '--monkhorst-pack',
                        metavar='K1,K2,K3',
                        help='Monkhorst-Pack sampling of BZ.  Example: ' +
                        '"4,4,4": 4x4x4 k-points, "4,4,4g": same set of ' +
                        'k-points shifted to include the Gamma point.')
        calc.add_option('--k-point-density', type='float', default=3.0,
                        help='Density of k-points in Angstrom.')
        calc.add_option('-p', '--parameters', metavar='key=value,...',
                        help='Comma-separated key=value pairs of ' +
                        'calculator specific parameters.')
        parser.add_option_group(calc)

    def parse(self, opts, args):
        mp = opts.monkhorst_pack
        if mp is not None:
            if mp[-1].lower() == 'g':
                kpts = np.array([int(k) for k in mp[:-1].split(',')])
                shift = 0.5 * ((kpts + 1) % 2) / kpts
                self.kpts = monkhorst_pack(kpts) + shift
            else:
                self.kpts = [int(k) for k in mp.split(',')]

        self.kptdensity = opts.k_point_density

        if opts.parameters:
            self.kwargs.update(str2dict(opts.parameters))


# Recognized names of calculators sorted alphabetically:
calcnames = ['abinit', 'aims', 'asap', 'castep', 'dftb', 'elk', 'emt',
             'exciting', 'fleur', 'gpaw', 'gaussian', 'hotbit', 'jacapo',
             'lammps', 'lj', 'mopac', 'morse',
             'nwchem', 'siesta', 'turbomole', 'vasp']

classnames = {'asap': 'EMT',
              'aims': 'Aims',
              'elk': 'ELK',
              'emt': 'EMT',
              'fleur': 'FLEUR',
              'gaussian': 'Gaussian',
              'jacapo': 'Jacapo',
              'lammps': 'LAMMPS',
              'lj': 'LennardJones',
              'mopac': 'Mopac',
              'morse': 'MorsePotential',
              'nwchem': 'NWChem',
              'vasp': 'Vasp'}


def calculator_factory(name, **kwargs):
    """Create an ASE calculator factory."""

    if name == 'abinit':
        from ase.calculators.abinit import Abinit
        return CalculatorFactory(Abinit, 'Abinit', 'label', **kwargs)

    if name == 'aims':
        from ase.calculators.aims import Aims
        return CalculatorFactory(Aims, 'aims', 'label', **kwargs)

    if name == 'nwchem':
        from ase.calculators.nwchem import NWChem
        return CalculatorFactory(NWChem, 'NWChem', 'label', 'no k-points',
                                 **kwargs)

    if name == 'asap':
        from asap3 import EMT
        return CalculatorFactory(EMT, 'Asap', None, 'no k-points', **kwargs)

    if name == 'elk':
        from ase.calculators.elk import ELK
        return CalculatorFactory(ELK, 'ELK', 'label', **kwargs)

    if name == 'fleur':
        from ase.calculators.fleur import FLEUR
        return CalculatorFactory(FLEUR, 'FLEUR', 'workdir', **kwargs)

    if name == 'gpaw':
        from gpaw.factory import GPAWFactory
        return GPAWFactory(**kwargs)

    if name == 'hotbit':
        from hotbit import Calculator
        return CalculatorFactory(Calculator, 'Hotbit', 'txt', 'no k-points',
                                 **kwargs)

    if name == 'jacapo':
        from ase.calculators.jacapo import Jacapo
        return CalculatorFactory(Jacapo, 'Jacapo', 'nc', **kwargs)

    if name == 'vasp':
        from ase.calculators.vasp import Vasp
        return CalculatorFactory(Vasp, 'Vasp', None, **kwargs)

    classname = classnames.get(name, name.title())
    module = __import__('ase.calculators.' + name, {}, None, [classname])
    Class = getattr(module, classname)

    if name in ['emt', 'gaussian', 'lammpsrun', 'lj', 'mopac', 'morse']:
        kpts = 'no k-points'
    else:
        kpts = None

    if name in ['emt', 'lj', 'morse']:
        label = None
    else:
        label = 'label'

    return CalculatorFactory(Class, classname, label, kpts, **kwargs)
