""" read and write gromacs geometry files
"""

from ase.atoms import Atoms
from ase.parallel import paropen


def read_gromacs(filename):
    """ From:
    http://manual.gromacs.org/current/online/gro.html
    C format
    "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" 
    python: starting from 0, including first excluding last
    0:4 5:10 10:15 15:20 20:28 28:36 36:44 44:52 52:60 60:68 

    Import gromacs geometry type files (.gro).
    Reads atom positions,
    velocities(if present) and
    simulation cell (if present)
    """

    from ase.data import chemical_symbols
    from ase import units

    atoms = Atoms()
    filed = open(filename, 'r')
    lines = filed.readlines()
    filed.close()
    positions = []
    gromacs_velocities = []
    symbols = []
    gromacs_residuenames = []
    gromacs_atomtypes = []
    for line in (lines[2:-1]):
        #print line[0:5]+':'+line[5:11]+':'+line[11:15]+':'+line[15:20]
        inp = line.split()
        floatvect = float(line[20:28]) * 10.0, \
            float(line[28:36]) * 10.0, \
            float(line[36:44]) * 10.0
        positions.append(floatvect)
        try:
            #velocities from nm/ps to ase units
            floatvect = \
                float(line[44:52]) * units.nm / (1000.0 * units.fs), \
                float(line[52:60]) * units.nm / (1000.0 * units.fs), \
                float(line[60:68]) * units.nm / (1000.0 * units.fs)
        except:
            floatvect = 0.0, 0.0, 0.0
        gromacs_velocities.append(floatvect)
        symbols.append(inp[1][0:2])
        gromacs_residuenames.append(inp[0])
        gromacs_atomtypes.append(inp[1])
    line = lines[-1]
    symbols_ok = []
    for isymbol in symbols:
        if isymbol in chemical_symbols:
            #ok atom name
            symbols_ok.append(isymbol)
        else:
            #not ok atom name
            symbols_ok.append(isymbol[0])
    atoms = Atoms(symbols_ok, positions)
    atoms.set_velocities(gromacs_velocities)

    if not atoms.has('residuenames'):
        atoms.new_array('residuenames', gromacs_residuenames, str)
        atoms.set_array('residuenames', gromacs_residuenames, str)
    if not atoms.has('atomtypes'):
        atoms.new_array('atomtypes', gromacs_atomtypes, str)
        atoms.set_array('atomtypes', gromacs_atomtypes, str)
    try:
        line = lines[-1]
        inp = line.split()
        floatvect0 = \
            float(inp[0]) * 10.0, \
            float(inp[1]) * 10.0, \
            float(inp[2]) * 10.0
        try:
            floatvect1 = \
                float(inp[3]) * 10.0, \
                float(inp[4]) * 10.0, \
                float(inp[5]) * 10.0
            floatvect2 = \
                float(inp[6]) * 10.0, \
                float(inp[7]) * 10.0, \
                float(inp[8]) * 10.0
            mycell = []
            #gromacs manual (manual.gromacs.org/online/gro.html) says:
            #v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
            #
            #v1(x) v2(y) v3(z) fv0[0 1 2]  v1(x) v2(x) v3(x)
            #v1(y) v1(z) v2(x) fv1[0 1 2]  v1(y) v2(y) v3(y)
            #v2(z) v3(x) v3(y) fv2[0 1 2]  v1(z) v2(z) v3(z)
            mycell += [[floatvect0[0], floatvect1[2], floatvect2[1]]]
            mycell += [[floatvect1[0], floatvect0[1], floatvect2[2]]]
            mycell += [[floatvect1[1], floatvect2[0], floatvect0[2]]]
            atoms.set_cell(mycell)
            atoms.set_pbc(True)
            print "9N"
        except:
            mycell = []
            #gromacs manual (manual.gromacs.org/online/gro.html) says:
            #v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
            mycell += [[floatvect0[0],           0.0,           0.0]]
            mycell += [[          0.0, floatvect0[1],           0.0]]
            mycell += [[          0.0,           0.0, floatvect0[2]]]
            atoms.set_cell(floatvect0)
            atoms.set_pbc(True)
    except:
        atoms.set_pbc(False)
    return atoms


def write_gromacs(fileobj, images):
    """Write gromacs geometry files (\*.gro).
    Writes:
    atom positions,
    velocities (if present, otherwise 0)
    and simulation cell (if present)
    """

    from ase import units

    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    natoms = images[-1].get_number_of_atoms()
    try:
        gromacs_residuenames = images[-1].get_array('residuenames')
    except:
        gromacs_residuenames = []
        for idum in range(natoms):
            gromacs_residuenames.append('1DUM')
    try:
        gromacs_atomtypes = images[-1].get_array('atomtypes')
    except:
        gromacs_atomtypes = images[-1].get_chemical_symbols()
    pos = images[-1].get_positions()
    pos = pos / 10.0
    try:
        vel = images[-1].get_velocities()
        vel = vel * 1000.0 * units.fs / units.nm
    except:
        vel = pos
        vel = pos * 0.0

    fileobj.write('#A Gromacs structure file written by ASE \n')
    fileobj.write('%5d \n' % images[-1].get_number_of_atoms())
    count = 1
    for resname, atomtype, xyz, vxyz in zip\
            (gromacs_residuenames, gromacs_atomtypes, pos, vel):
        fileobj.write(\
            '   %5s  %5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n' % \
                (resname, atomtype, count, \
                xyz[0], xyz[1], xyz[2], \
                vxyz[0], vxyz[1], vxyz[2]))
        count = count + 1
    if images[-1].get_pbc().any():
        mycell = images[-1].get_cell()
        #gromacs manual (manual.gromacs.org/online/gro.html) says:
        #v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
        #
        #cell[0,0] cell[1,0] cell[2,0] v1(x) v2(y) v3(z) fv0[0 1 2]
        #cell[0,1] cell[1,1] cell[2,1] v1(y) v1(z) v2(x) fv1[0 1 2]
        #cell[0,2] cell[1,2] cell[2,2] v2(z) v3(x) v3(y) fv2[0 1 2]
        fileobj.write('%10.5f%10.5f%10.5f' \
                          % (mycell[0, 0] * 0.1, \
                                 mycell[1, 1] * 0.1, \
                                 mycell[2, 2] * 0.1))
        fileobj.write('%10.5f%10.5f%10.5f' \
                          % (mycell[1, 0] * 0.1, \
                                 mycell[2, 0] * 0.1, \
                                 mycell[0, 1] * 0.1))
        fileobj.write('%10.5f%10.5f%10.5f\n' \
                          % (mycell[2, 1] * 0.1, \
                                 mycell[0, 2] * 0.1, \
                                 mycell[1, 2] * 0.1))
    return
