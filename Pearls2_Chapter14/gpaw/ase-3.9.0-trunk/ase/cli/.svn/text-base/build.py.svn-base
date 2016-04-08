import sys
import optparse
    
import numpy as np

from ase.db import connect
from ase.lattice import bulk
from ase.io import read, write
from ase.visualize import view
from ase.structure import molecule
from ase.atoms import Atoms, string2symbols
from ase.data import ground_state_magnetic_moments
from ase.data import atomic_numbers, covalent_radii


def main():
    parser = optparse.OptionParser(
        usage='%prog [options] name/input-file [output-file]')
    add = parser.add_option
    add('-M', '--magnetic-moment',
        metavar='M1,M2,...',
        help='Magnetic moment(s).  ' +
        'Use "-M 1" or "-M 2.3,-2.3".')
    add('--modify', metavar='...',
        help='Modify atoms with Python statement.  ' +
        'Example: --modify="atoms.positions[-1,2]+=0.1".')
    add('-v', '--vacuum', type=float, default=3.0,
        help='Amount of vacuum to add around isolated atoms '
        '(in Angstrom).')
    add('--unit-cell',
        help='Unit cell.  Examples: "10.0" or "9,10,11" ' +
        '(in Angstrom).')
    add('--bond-length', type=float,
        help='Bond length of dimer in Angstrom.')
    add('-x', '--crystal-structure',
        help='Crystal structure.',
        choices=['sc', 'fcc', 'bcc', 'hcp', 'diamond',
                 'zincblende', 'rocksalt', 'cesiumchloride',
                 'fluorite', 'wurtzite'])
    add('-a', '--lattice-constant', default='',
        help='Lattice constant(s) in Angstrom.')
    add('--orthorhombic', action='store_true',
        help='Use orthorhombic unit cell.')
    add('--cubic', action='store_true',
        help='Use cubic unit cell.')
    add('-r', '--repeat',
        help='Repeat unit cell.  Use "-r 2" or "-r 2,3,1".')
    add('-g', '--gui', action='store_true')

    opts, args = parser.parse_args()
    if len(args) == 0 or len(args) > 2:
        parser.error('Wrong number of arguments!')

    name = args.pop(0)

    if '.' in name:
        # Read from file:
        atoms = read(name)
    elif opts.crystal_structure:
        atoms = build_bulk(name, opts)
    else:
        atoms = build_molecule(name, opts)

    if opts.magnetic_moment:
        magmoms = np.array(
            [float(m) for m in opts.magnetic_moment.split(',')])
        atoms.set_initial_magnetic_moments(
            np.tile(magmoms, len(atoms) // len(magmoms)))

    if opts.modify:
        exec opts.modify in {'atoms': atoms}

    if opts.repeat is not None:
        r = opts.repeat.split(',')
        if len(r) == 1:
            r = 3 * r
        atoms = atoms.repeat([int(c) for c in r])

    if opts.gui:
        view(atoms)
        
    if args:
        write(args[0], atoms)
    elif sys.stdout.isatty():
        write(name + '.json', atoms)
    else:
        con = connect(sys.stdout, type='json')
        con.write(atoms, name=name)


def build_molecule(name, opts):
    try:
        # Known molecule or atom?
        atoms = molecule(name)
    except NotImplementedError:
        symbols = string2symbols(name)
        if len(symbols) == 1:
            Z = atomic_numbers[symbols[0]]
            magmom = ground_state_magnetic_moments[Z]
            atoms = Atoms(name, magmoms=[magmom])
        elif len(symbols) == 2:
            # Dimer
            if opts.bond_length is None:
                b = (covalent_radii[atomic_numbers[symbols[0]]] +
                     covalent_radii[atomic_numbers[symbols[1]]])
            else:
                b = opts.bond_length
            atoms = Atoms(name, positions=[(0, 0, 0),
                                           (b, 0, 0)])
        else:
            raise ValueError('Unknown molecule: ' + name)
    else:
        if len(atoms) == 2 and opts.bond_length is not None:
            atoms.set_distance(0, 1, opts.bond_length)

    if opts.unit_cell is None:
        atoms.center(vacuum=opts.vacuum)
    else:
        a = [float(x) for x in opts.unit_cell.split(',')]
        if len(a) == 1:
            cell = [a[0], a[0], a[0]]
        elif len(a) == 3:
            cell = a
        else:
            a, b, c, alpha, beta, gamma = a
            degree = np.pi / 180.0
            cosa = np.cos(alpha * degree)
            cosb = np.cos(beta * degree)
            sinb = np.sin(beta * degree)
            cosg = np.cos(gamma * degree)
            sing = np.sin(gamma * degree)
            cell = [[a, 0, 0],
                    [b * cosg, b * sing, 0],
                    [c * cosb, c * (cosa - cosb * cosg) / sing,
                     c * np.sqrt(
                        sinb**2 - ((cosa - cosb * cosg) / sing)**2)]]
        atoms.cell = cell
        atoms.center()

    return atoms


def build_bulk(name, opts):
    L = opts.lattice_constant.replace(',', ' ').split()
    d = dict([(key, float(x)) for key, x in zip('ac', L)])
    atoms = bulk(name, crystalstructure=opts.crystal_structure,
                 a=d.get('a'), c=d.get('c'),
                 orthorhombic=opts.orthorhombic, cubic=opts.cubic)

    M, X = {'Fe': (2.3, 'bcc'),
            'Co': (1.2, 'hcp'),
            'Ni': (0.6, 'fcc')}.get(name, (None, None))
    if M is not None and opts.crystal_structure == X:
        atoms.set_initial_magnetic_moments([M] * len(atoms))

    return atoms
