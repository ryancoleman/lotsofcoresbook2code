"""Reads quantum espresso files. Tested for output on PWSCF v.5.0.2, only
for typical output of input files made with ASE -- that is, ibrav=0."""

import numpy as np
from ase.atoms import Atoms, Atom
from ase import units
from ase.calculators.singlepoint import SinglePointCalculator


def read_espresso_out(fileobj, index):
    """Reads quantum espresso output text files."""
    if isinstance(fileobj, str):
        fileobj = open(fileobj, 'rU')
    lines = fileobj.readlines()
    images = []

    # Get unit cell info.
    bl_line = [line for line in lines if 'bravais-lattice index' in line]
    if len(bl_line) != 1:
        raise NotImplementedError('Unsupported: unit cell changing.')
    bl_line = bl_line[0].strip()
    brav_latt_index = bl_line.split('=')[1].strip()
    if brav_latt_index != '0':
        raise NotImplementedError('Supported only for Bravais-lattice '
                                  'index of 0 (free).')
    lp_line = [line for line in lines if 'lattice parameter (alat)' in
               line]
    if len(lp_line) != 1:
        raise NotImplementedError('Unsupported: unit cell changing.')
    lp_line = lp_line[0].strip().split('=')[1].strip().split()[0]
    lattice_parameter = float(lp_line) * units.Bohr
    ca_line_no = [number for (number, line) in enumerate(lines) if
                  'crystal axes: (cart. coord. in units of alat)' in line]
    if len(ca_line_no) != 1:
        raise NotImplementedError('Unsupported: unit cell changing.')
    ca_line_no = int(ca_line_no[0])
    cell = np.zeros((3, 3))
    for number, line in enumerate(lines[ca_line_no + 1: ca_line_no + 4]):
        line = line.split('=')[1].strip()[1:-1]
        values = [eval(value) for value in line.split()]
        cell[number, 0] = values[0]
        cell[number, 1] = values[1]
        cell[number, 2] = values[2]
    cell *= lattice_parameter

    # Find atomic positions and add to images.
    for number, line in enumerate(lines):
        key = 'Begin final coordinates'  # these just reprint last posn.
        if key in line:
            break
        key = 'Cartesian axes'
        if key in line:
            atoms = make_atoms(number, lines, key, cell)
            images.append(atoms)
        key = 'ATOMIC_POSITIONS (crystal)'
        if key in line:
            atoms = make_atoms(number, lines, key, cell)
            images.append(atoms)
    return images[index]


def make_atoms(index, lines, key, cell):
    """Scan through lines to get the atomic positions."""
    atoms = Atoms()
    if key == 'Cartesian axes':
        for line in lines[index + 3:]:
            entries = line.split()
            if len(entries) == 0:
                break
            symbol = entries[1][:-1]
            x = float(entries[6])
            y = float(entries[7])
            z = float(entries[8])
            atoms.append(Atom(symbol, (x, y, z)))
        atoms.set_cell(cell)
    elif key == 'ATOMIC_POSITIONS (crystal)':
        for line in lines[index + 1:]:
            entries = line.split()
            if len(entries) == 0 or (entries[0] == 'End'):
                break
            symbol = entries[0][:-1]
            x = float(entries[1])
            y = float(entries[2])
            z = float(entries[3])
            atoms.append(Atom(symbol, (x, y, z)))
        atoms.set_cell(cell, scale_atoms=True)
    # Energy is located after positions.
    energylines = [number for number, line in enumerate(lines) if
                   ('!' in line and 'total energy' in line)]
    energyline = min([n for n in energylines if n > index])
    energy = float(lines[energyline].split()[-2]) * units.Ry
    # Forces are located after positions.
    forces = np.zeros((len(atoms), 3))
    forcelines = [number for number, line in enumerate(lines) if
                  'Forces acting on atoms (Ry/au):' in line]
    forceline = min([n for n in forcelines if n > index])
    for line in lines[forceline + 4:]:
        words = line.split()
        if len(words) == 0:
            break
        fx = float(words[-3])
        fy = float(words[-2])
        fz = float(words[-1])
        atom_number = int(words[1]) - 1
        forces[atom_number] = (fx, fy, fz)
    forces *= units.Ry / units.Bohr
    calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
    atoms.set_calculator(calc)
    return atoms


def read_espresso_in(fileobj):
    """Reads espresso input files."""
    if isinstance(fileobj, str):
        fileobj = open(fileobj, 'rU')
    data, extralines = read_fortran_namelist(fileobj)
    positions, method = get_atomic_positions(extralines,
                                             n_atoms=data['system']['nat'])
    cell = get_cell_parameters(extralines)
    if data['system']['ibrav'] == 0:
        atoms = build_atoms(positions, method, cell,
                            data['system']['celldm(1)'])
    else:
        raise NotImplementedError('ibrav=%i not implemented.' %
                                  data['system']['ibrav'])
    return atoms


def build_atoms(positions, method, cell, alat):
    """Creates the atoms for a quantum espresso in file."""
    if method != 'crystal':
        raise NotImplementedError('Only supported for crystal method of '
                                  'ATOMIC_POSITIONS, not %s.' % method)
    atoms = Atoms()
    for el, (x, y, z) in positions:
        atoms.append(Atom(el, (x, y, z)))
    cell *= f2f(alat) * units.Bohr
    atoms.set_cell(cell, scale_atoms=True)
    return atoms


def get_atomic_positions(lines, n_atoms):
    """Returns the atomic positions of the atoms as an (ordered) list from
    the lines of text of the espresso input file."""
    atomic_positions = []
    line = [n for (n, l) in enumerate(lines) if 'ATOMIC_POSITIONS' in l]
    if len(line) == 0:
        return None
    if len(line) > 1:
        raise RuntimeError('More than one ATOMIC_POSITIONS section?')
    line_no = line[0]
    for line in lines[line_no + 1:line_no + n_atoms + 1]:
        el, x, y, z = line.split()
        atomic_positions.append([el[:-1], (f2f(x), f2f(y), f2f(z))])
    line = lines[line_no]
    if '{' in line:
        method = line[line.find('{') + 1:line.find('}')]
    elif '(' in line:
        method = line[line.find('(') + 1:line.find(')')]
    else:
        method = None
    return atomic_positions, method


def get_cell_parameters(lines):
    """Returns the cell parameters as a matrix."""
    cell_parameters = np.zeros((3, 3))
    line = [n for (n, l) in enumerate(lines) if 'CELL_PARAMETERS' in l]
    if len(line) == 0:
        return None
    if len(line) > 1:
        raise RuntimeError('More than one CELL_PARAMETERS section?')
    line_no = line[0]
    for vector, line in enumerate(lines[line_no + 1:line_no + 4]):
        x, y, z = line.split()
        cell_parameters[vector] = (f2f(x), f2f(y), f2f(z))
    return cell_parameters


def read_fortran_namelist(fileobj):
    """Takes a fortran-namelist formatted file and returns appropriate
    dictionaries, followed by lines of text that do not fit this pattern.
    """
    data = {}
    extralines = []
    indict = False
    fileobj.seek(0)
    for line in fileobj.readlines():
        if indict and line.strip().startswith('/'):
            indict = False
        elif line.strip().startswith('&'):
            indict = True
            dictname = line.strip()[1:].lower()
            data[dictname] = {}
        elif (not indict) and (len(line.strip()) > 0):
            extralines.append(line)
        elif indict:
            key, value = line.strip().split('=')
            if value.endswith(','):
                value = value[:-1]
            value = value.strip()
            try:
                value = eval(value)
            except SyntaxError:
                value = {'.true.': True, '.false.': False}.get(value, value)
            data[dictname][key.strip()] = value
    return data, extralines


def f2f(value):
    """Converts a fortran-formatted double precision number (e.g., 2.323d2)
    to a python float. value should be a string."""
    value = value.replace('d', 'e')
    value = value.replace('D', 'e')
    return float(value)
