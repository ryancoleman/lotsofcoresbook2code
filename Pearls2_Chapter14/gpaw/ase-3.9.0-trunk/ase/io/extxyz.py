"""
Extended XYZ support

Read/write files in "extended" XYZ format, storing additional
per-configuration information as key-value pairs on the XYZ
comment line, and additional per-atom properties as extra columns.

See http://jrkermode.co.uk/quippy/io.html#extendedxyz for a full
description of the Extended XYZ file format.

Contributed by James Kermode <james.kermode@gmail.com>
"""

import re
import numpy as np

from ase.atoms import Atoms
from ase.parallel import paropen
from ase.calculators.calculator import all_properties
from ase.calculators.singlepoint import SinglePointCalculator

__all__ = ['read_xyz', 'write_xyz']

PROPERTY_NAME_MAP = {'positions': 'pos',
                     'numbers': 'Z',
                     'charges': 'charge',
                     'symbols': 'species'}

REV_PROPERTY_NAME_MAP = dict(zip(PROPERTY_NAME_MAP.values(),
                                 PROPERTY_NAME_MAP.keys()))

KEY_QUOTED_VALUE = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)' +
                              r'\s*=\s*["\{\}]([^"\{\}]+)["\{\}e+-]\s*')
KEY_VALUE = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*=' +
                       r'\s*([-0-9A-Za-z_.:\[\]()e+-]+)\s*')
KEY_RE = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*')


def key_val_str_to_dict(s):
    """
    Parse strings in the form 'key1=value1 key2="quoted value"'
    """
    d = {}
    s = s.strip()
    while True:
        # Match quoted string first, then fall through to plain key=value
        m = KEY_QUOTED_VALUE.match(s)
        if m is None:
            m = KEY_VALUE.match(s)
            if m is not None:
                s = KEY_VALUE.sub('', s, 1)
            else:
                # Just a key with no value
                m = KEY_RE.match(s)
                if m is not None:
                    s = KEY_RE.sub('', s, 1)
        else:
            s = KEY_QUOTED_VALUE.sub('', s, 1)

        if m is None:
            break        # No more matches

        key = m.group(1)
        try:
            value = m.group(2)
        except IndexError:
            # default value is 'T' (True)
            value = 'T'

        # Try to convert to (arrays of) floats, ints
        try:
            numvalue = []
            for x in value.split():
                if x.find('.') == -1:
                    numvalue.append(int(float(x)))
                else:
                    numvalue.append(float(x))
            if len(numvalue) == 1:
                numvalue = numvalue[0]         # Only one number
            elif len(numvalue) == 9:
                # special case: 3x3 matrix, fortran ordering
                numvalue = np.array(numvalue).reshape((3, 3), order='F')
            else:
                numvalue = np.array(numvalue)  # vector
            value = numvalue
        except ValueError:
            pass

        # Parse boolean values: 'T' -> True, 'F' -> False,
        #                       'T T F' -> [True, True, False]
        if isinstance(value, str):
            str_to_bool = {'T': True, 'F': False}

            if len(value.split()) > 1:
                if all([x in str_to_bool.keys() for x in value.split()]):
                    value = [str_to_bool[x] for x in value.split()]
            elif value in str_to_bool:
                value = str_to_bool[value]

        d[key] = value

    return d


def key_val_dict_to_str(d, sep=' '):
    """
    Convert atoms.info dictionary to extended XYZ string representation
    """
    if len(d) == 0:
        return ''
    s = ''
    type_val_map = {(bool, True): 'T',
                    (bool, False): 'F',
                    (np.bool_, True): 'T',
                    (np.bool_, False): 'F'}

    s = ''
    for key in d.keys():
        val = d[key]

        if hasattr(val, '__iter__'):
            val = np.array(val)
            val = ' '.join(str(type_val_map.get((type(x), x), x))
                           for x in val.reshape(val.size, order='F'))
            val.replace('[', '')
            val.replace(']', '')
        else:
            val = type_val_map.get((type(val), val), val)

        if val is None:
            s = s + '%s%s' % (key, sep)
        elif isinstance(val, basestring) and ' ' in val:
            s = s + '%s="%s"%s' % (key, val, sep)
        else:
            s = s + '%s=%s%s' % (key, str(val), sep)

    return s.strip()


def parse_properties(prop_str):
    """
    Parse extended XYZ properties format string

    Format is "[NAME:TYPE:NCOLS]...]", e.g. "species:S:1:pos:R:3".
    NAME is the name of the property.
    TYPE is one of R, I, S, L for real, integer, string and logical.
    NCOLS is number of columns for that property.
    """

    properties = {}
    properties_list = []
    dtypes = []
    converters = []

    fields = prop_str.split(':')

    def parse_bool(x):
        """
        Parse bool to string
        """
        return {'T': True, 'F': False,
                'True': True, 'False': False}.get(x)

    fmt_map = {'R': ('d', float),
               'I': ('i', int),
               'S': (object, str),
               'L': ('bool', parse_bool)}

    for name, ptype, cols in zip(fields[::3],
                                 fields[1::3],
                                 [int(x) for x in fields[2::3]]):
        if ptype not in ('R', 'I', 'S', 'L'):
            raise ValueError('Unknown property type: ' + ptype)
        ase_name = REV_PROPERTY_NAME_MAP.get(name, name)

        dtype, converter = fmt_map[ptype]
        if cols == 1:
            dtypes.append((name, dtype))
            converters.append(converter)
        else:
            for c in range(cols):
                dtypes.append((name + str(c), dtype))
                converters.append(converter)

        properties[name] = (ase_name, cols)
        properties_list.append(name)

    dtype = np.dtype(dtypes)
    return properties, properties_list, dtype, converters


def read_xyz(fileobj, index=-1):
    """
    Read from a file in Extended XYZ format

    index is the frame to read, default is last frame (index=-1).
    """
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    if not isinstance(index, int) and not isinstance(index, slice):
        raise TypeError('Index argument is neither slice nor integer!')

    fileobj.seek(0)
    natoms = int(fileobj.readline())

    for i, l in enumerate(fileobj):
        pass

    ln = i + 2
    lnsnp = natoms + 2
    lastsnap = ln // lnsnp

    rvrs = False

    if isinstance(index, int):
        if index < 0:
            tmpsnp = lastsnap + index
            trbl = range(tmpsnp, tmpsnp + 1, 1)
        else:
            trbl = range(index, index + 1, 1)
        rtnndx = -1
    elif isinstance(index, slice):
        start = index.start
        stop = index.stop
        step = index.step

        if start is None:
            start = 0
        elif start < 0:
            start = lastsnap + start

        if step is None:
            step = 1

        if stop is None:
            stop = lastsnap
        elif stop < 0:
            stop = lastsnap + stop

        trbl = range(start, stop, step)

        if step < 0:
            rvrs = True
            trbl.reverse()

        rtnndx = slice(len(trbl))

    images = []
    current = 0

    fileobj.seek(0)

    for index in trbl:
        for lnndx in range(current, index * lnsnp):
            line = fileobj.readline()

        line = fileobj.readline()
        line = fileobj.readline()

        info = key_val_str_to_dict(line)

        pbc = None
        if 'pbc' in info:
            pbc = info['pbc']
            del info['pbc']
        elif 'Lattice' in info:
            # default pbc for extxyz file containing Lattice
            # is True in all directions
            pbc = [True, True, True]

        cell = None
        if 'Lattice' in info:
            # NB: ASE cell is transpose of extended XYZ lattice
            cell = info['Lattice'].T
            del info['Lattice']

        if 'Properties' not in info:
            # Default set of properties is atomic symbols and positions only
            info['Properties'] = 'species:S:1:pos:R:3'
        properties, names, dtype, convs = parse_properties(info['Properties'])
        del info['Properties']

        data = []
        for ln in range(natoms):
            line = fileobj.readline()
            vals = line.split()
            row = tuple([conv(val) for conv, val in zip(convs, vals)])
            data.append(row)

        try:
            data = np.array(data, dtype)
        except TypeError:
            raise IOError('Badly formatted data, ' +
                          'or end of file reached before end of frame')

        arrays = {}
        for name in names:
            ase_name, cols = properties[name]
            if cols == 1:
                value = data[name]
            else:
                value = np.vstack([data[name + str(c)]
                                   for c in range(cols)]).T
            arrays[ase_name] = value

        symbols = None
        if 'symbols' in arrays:
            symbols = arrays['symbols']
            del arrays['symbols']

        numbers = None
        duplicate_numbers = None
        if 'numbers' in arrays:
            if symbols is None:
                numbers = arrays['numbers']
            else:
                duplicate_numbers = arrays['numbers']
            del arrays['numbers']

        positions = None
        if 'positions' in arrays:
            positions = arrays['positions']
            del arrays['positions']

        atoms = Atoms(symbols=symbols,
                      positions=positions,
                      numbers=numbers,
                      cell=cell,
                      pbc=pbc,
                      info=info)

        for (name, array) in arrays.iteritems():
            atoms.new_array(name, array)

        if duplicate_numbers is not None:
            atoms.set_atomic_numbers(duplicate_numbers)

        # Load results of previous calculations into SinglePointCalculator
        results = {}
        for key in atoms.info.keys():
            if key in all_properties:
                results[key] = atoms.info[key]
        for key in atoms.arrays.keys():
            if key in all_properties:
                results[key] = atoms.arrays[key]
        if results != {}:
            calculator = SinglePointCalculator(atoms, **results)
            atoms.set_calculator(calculator)

        images.append(atoms)
        current = (index + 1) * lnsnp

    if rvrs:
        images.reverse()

    return images[rtnndx]


def output_column_format(atoms, columns, arrays, write_info=True):
    """
    Helper function to build extended XYZ comment line
    """
    fmt_map = {'d': ('R', '%16.8f '),
               'f': ('R', '%16.8f '),
               'i': ('I', '%8d '),
               'O': ('S', '%s'),
               'S': ('S', '%s'),
               'b': ('L', ' %.1s ')}

    # NB: Lattice is stored as tranpose of ASE cell,
    # with Fortran array ordering
    lattice_str = ('Lattice="' +
                   ' '.join([str(x) for x in np.reshape(atoms.cell.T,
                                                        9, order='F')]) +
                   '"')

    property_names = []
    property_types = []
    property_ncols = []
    dtypes = []
    formats = []

    for column in columns:
        array = arrays[column]
        dtype = array.dtype

        property_name = PROPERTY_NAME_MAP.get(column, column)
        property_type, fmt = fmt_map[dtype.kind]
        property_names.append(property_name)
        property_types.append(property_type)

        if len(array.shape) == 1:
            ncol = 1
            dtypes.append((column, dtype))
        else:
            ncol = array.shape[1]
            for c in range(ncol):
                dtypes.append((column + str(c), dtype))

        formats.extend([fmt] * ncol)
        property_ncols.append(ncol)

    props_str = ':'.join([':'.join(x) for x in
                          zip(property_names,
                              property_types,
                              [str(nc) for nc in property_ncols])])

    comment = lattice_str + ' Properties=' + props_str
    info = {}
    if write_info:
        info.update(atoms.info)
    info['pbc'] = atoms.get_pbc()  # always save periodic boundary conditions
    comment += ' ' + key_val_dict_to_str(info)

    dtype = np.dtype(dtypes)
    fmt = ''.join(formats) + '\n'

    return comment, property_ncols, dtype, fmt


def write_xyz(fileobj, images, columns=None, write_info=True):
    """
    Write output in extended XYZ format

    Optionally, specify which columns (arrays) to include in output,
    and whether to write the contents of the Atoms.info dict to the
    XYZ comment line (default is True)
    """
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    for atoms in images:
        natoms = len(atoms)

        if columns is None:
            columns = ['symbols'] + atoms.arrays.keys()

        # Move symbols and positions to first two properties
        if 'symbols' in columns:
            i = columns.index('symbols')
            columns[0], columns[i] = columns[i], columns[0]

        if 'positions' in columns:
            i = columns.index('positions')
            columns[1], columns[i] = columns[i], columns[1]

        # Check first column "looks like" atomic symbols
        if columns[0] in atoms.arrays:
            symbols = atoms.arrays[columns[0]]
        else:
            symbols = atoms.get_chemical_symbols()
        if not isinstance(symbols[0], basestring):
            raise ValueError('First column must be symbols-like')

        # Check second column "looks like" atomic positions
        pos = atoms.arrays[columns[1]]
        if pos.shape != (natoms, 3) or pos.dtype.kind != 'f':
            raise ValueError('Second column must be position-like')

        # Collect data to be written out
        arrays = {}
        for column in columns:
            if column in atoms.arrays:
                arrays[column] = atoms.arrays[column]
            elif column == 'symbols':
                arrays[column] = np.array(symbols)
            else:
                raise ValueError('Missing array %s')

        comm, ncols, dtype, fmt = output_column_format(atoms, columns, arrays,
                                                       write_info)

        # Pack columns into record array
        data = np.zeros(natoms, dtype)
        for column, ncol in zip(columns, ncols):
            value = arrays[column]
            if ncol == 1:
                data[column] = value
            else:
                for c in range(ncol):
                    data[column + str(c)] = value[:, c]

        # Write the output
        fileobj.write('%d\n' % natoms)
        fileobj.write('%s\n' % comm)
        for i in range(natoms):
            fileobj.write(fmt % tuple(data[i]))
