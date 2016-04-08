import os
import sys
from tarfile import is_tarfile
from zipfile import is_zipfile

from ase.atoms import Atoms
from ase.units import Bohr, Hartree
from ase.io.trajectory import PickleTrajectory
from ase.io.bundletrajectory import BundleTrajectory
from ase.io.netcdftrajectory import NetCDFTrajectory
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.calculators.singlepoint import SinglePointKPoint

__all__ = ['read', 'write', 'PickleTrajectory', 'BundleTrajectory',
           'NetCDFTrajectory']


def read(filename, index=None, format=None):
    """Read Atoms object(s) from file.

    filename: str
        Name of the file to read from.
    index: int or slice
        If the file contains several configurations, the last configuration
        will be returned by default.  Use index=n to get configuration
        number n (counting from zero).
    format: str
        Used to specify the file-format.  If not given, the
        file-format will be guessed by the *filetype* function.

    Known formats:

    =========================  =============
    format                     short name
    =========================  =============
    GPAW restart-file          gpw
    Dacapo netCDF output file  dacapo
    Old ASE netCDF trajectory  nc
    Virtual Nano Lab file      vnl
    ASE pickle trajectory      traj
    ASE bundle trajectory      bundle
    GPAW text output           gpaw-text
    CUBE file                  cube
    XCrySDen Structure File    xsf
    Dacapo text output         dacapo-text
    XYZ-file                   xyz
    VASP POSCAR/CONTCAR file   vasp
    VASP OUTCAR file           vasp_out
    VASP XDATCAR file          vasp_xdatcar
    SIESTA STRUCT file         struct_out
    ABINIT input file          abinit
    V_Sim ascii file           v_sim
    Protein Data Bank          pdb
    CIF-file                   cif
    FHI-aims geometry file     aims
    FHI-aims output file       aims_out
    VTK XML Image Data         vti
    VTK XML Structured Grid    vts
    VTK XML Unstructured Grid  vtu
    TURBOMOLE coord file       tmol
    TURBOMOLE gradient file    tmol-gradient
    exciting input             exi
    AtomEye configuration      cfg
    WIEN2k structure file      struct
    DftbPlus input file        dftb
    CASTEP geom file           cell
    CASTEP output file         castep
    CASTEP trajectory file     geom
    ETSF format                etsf.nc
    DFTBPlus GEN format        gen
    CMR db/cmr-file            db
    CMR db/cmr-file            cmr
    LAMMPS dump file           lammps
    EON reactant.con file      eon
    Gromacs coordinates        gro
    Gaussian com (input) file  gaussian
    Gaussian output file       gaussian_out
    Quantum espresso in file   esp_in
    Quantum espresso out file  esp_out
    Extended XYZ file          extxyz
    NWChem input file          nw
    Materials Studio file      xsd
    =========================  =============

    """
    if isinstance(filename, str) and (
        '.json@' in filename or
        '.db@' in filename or
        filename.startswith('pg://') and '@' in filename):
        filename, index = filename.rsplit('@', 1)
        if index.isdigit():
            index = int(index)
    else:
        if isinstance(filename, str):
            p = filename.rfind('@')
            if p != -1:
                try:
                    index = string2index(filename[p + 1:])
                except ValueError:
                    pass
                else:
                    filename = filename[:p]

        if isinstance(index, str):
            index = string2index(index)

    if format is None:
        format = filetype(filename)

    if format.startswith('gpw'):
        import gpaw
        r = gpaw.io.open(filename, 'r')
        positions = r.get('CartesianPositions') * Bohr
        numbers = r.get('AtomicNumbers')
        cell = r.get('UnitCell') * Bohr
        pbc = r.get('BoundaryConditions')
        tags = r.get('Tags')
        magmoms = r.get('MagneticMoments')
        energy = r.get('PotentialEnergy') * Hartree

        if r.has_array('CartesianForces'):
            forces = r.get('CartesianForces') * Hartree / Bohr
        else:
            forces = None

        atoms = Atoms(positions=positions,
                      numbers=numbers,
                      cell=cell,
                      pbc=pbc)
        if tags.any():
            atoms.set_tags(tags)

        if magmoms.any():
            atoms.set_initial_magnetic_moments(magmoms)
        else:
            magmoms = None

        atoms.calc = SinglePointDFTCalculator(atoms, energy=energy,
                                              forces=forces, magmoms=magmoms)
        kpts = []
        if r.has_array('IBZKPoints'):
            for w, kpt, eps_n, f_n in zip(r.get('IBZKPointWeights'),
                                          r.get('IBZKPoints'),
                                          r.get('Eigenvalues'),
                                          r.get('OccupationNumbers')):
                kpts.append(SinglePointKPoint(w, kpt[0], kpt[1],
                                              eps_n[0], f_n[0]))
        atoms.calc.kpts = kpts

        return atoms

    if format in ['json', 'db', 'postgresql']:
        from ase.db.core import connect, dict2atoms
        if index == slice(None, None):
            index = None
        images = [dict2atoms(d)
                  for d in connect(filename, format).select(index)]
        if len(images) == 1:
            return images[0]
        else:
            return images

    if index is None:
        index = -1
        
    if format == 'castep':
        from ase.io.castep import read_castep
        return read_castep(filename, index)

    if format == 'castep_cell':
        import ase.io.castep
        return ase.io.castep.read_cell(filename, index)

    if format == 'castep_geom':
        import ase.io.castep
        return ase.io.castep.read_geom(filename, index)

    if format == 'exi':
        from ase.io.exciting import read_exciting
        return read_exciting(filename, index)

    if format in ['xyz', 'extxyz']:
        from ase.io.extxyz import read_xyz
        return read_xyz(filename, index)

    if format == 'traj':
        from ase.io.trajectory import read_trajectory
        return read_trajectory(filename, index)

    if format == 'bundle':
        from ase.io.bundletrajectory import read_bundletrajectory
        return read_bundletrajectory(filename, index)

    if format == 'cube':
        from ase.io.cube import read_cube
        return read_cube(filename, index)

    if format == 'nc':
        from ase.io.netcdf import read_netcdf
        return read_netcdf(filename, index)

    if format == 'gpaw-text':
        from ase.io.gpawtext import read_gpaw_text
        return read_gpaw_text(filename, index)

    if format == 'dacapo-text':
        from ase.io.dacapo import read_dacapo_text
        return read_dacapo_text(filename)

    if format == 'dacapo':
        from ase.io.dacapo import read_dacapo
        return read_dacapo(filename)

    if format == 'xsf':
        from ase.io.xsf import read_xsf
        return read_xsf(filename, index)

    if format == 'vasp':
        from ase.io.vasp import read_vasp
        return read_vasp(filename)

    if format == 'vasp_out':
        from ase.io.vasp import read_vasp_out
        return read_vasp_out(filename, index)

    if format == 'vasp_xdatcar':
        from ase.io.vasp import read_vasp_xdatcar
        return read_vasp_xdatcar(filename, index)

    if format == 'abinit':
        from ase.io.abinit import read_abinit
        return read_abinit(filename)

    if format == 'v_sim':
        from ase.io.v_sim import read_v_sim
        return read_v_sim(filename)

    if format == 'mol':
        from ase.io.mol import read_mol
        return read_mol(filename)

    if format == 'pdb':
        from ase.io.pdb import read_pdb
        return read_pdb(filename, index)

    if format == 'cif':
        from ase.io.cif import read_cif
        return read_cif(filename, index)

    if format == 'struct':
        from ase.io.wien2k import read_struct
        return read_struct(filename)

    if format == 'struct_out':
        from ase.io.siesta import read_struct
        return read_struct(filename)

    if format == 'vti':
        from ase.io.vtkxml import read_vti
        return read_vti(filename)

    if format == 'vts':
        from ase.io.vtkxml import read_vts
        return read_vts(filename)

    if format == 'vtu':
        from ase.io.vtkxml import read_vtu
        return read_vtu(filename)

    if format == 'aims':
        from ase.io.aims import read_aims
        return read_aims(filename)

    if format == 'aims_out':
        from ase.io.aims import read_aims_output
        return read_aims_output(filename, index)

    if format == 'iwm':
        from ase.io.iwm import read_iwm
        return read_iwm(filename)

    if format == 'Cmdft':
        from ase.io.cmdft import read_I_info
        return read_I_info(filename)

    if format == 'tmol':
        from ase.io.turbomole import read_turbomole
        return read_turbomole(filename)

    if format == 'tmol-gradient':
        from ase.io.turbomole import read_turbomole_gradient
        return read_turbomole_gradient(filename)

    if format == 'cfg':
        from ase.io.cfg import read_cfg
        return read_cfg(filename)

    if format == 'dftb':
        from ase.io.dftb import read_dftb
        return read_dftb(filename)

    if format == 'sdf':
        from ase.io.sdf import read_sdf
        return read_sdf(filename)

    if format == 'etsf':
        from ase.io.etsf import ETSFReader
        return ETSFReader(filename).read_atoms()

    if format == 'gen':
        from ase.io.gen import read_gen
        return read_gen(filename)

    if format == 'cmr':
        from ase.io.cmr_io import read_db
        return read_db(filename, index)

    if format == 'lammps':
        from ase.io.lammpsrun import read_lammps_dump
        return read_lammps_dump(filename, index)

    if format == 'eon':
        from ase.io.eon import read_reactant_con
        return read_reactant_con(filename)

    if format == 'gromacs':
        from ase.io.gromacs import read_gromacs
        return read_gromacs(filename)

    if format == 'gaussian':
        from ase.io.gaussian import read_gaussian
        return read_gaussian(filename)

    if format == 'gaussian_out':
        from ase.io.gaussian import read_gaussian_out
        return read_gaussian_out(filename, index)

    if format == 'esp_in':
        from ase.io.espresso import read_espresso_in
        return read_espresso_in(filename)

    if format == 'esp_out':
        from ase.io.espresso import read_espresso_out
        return read_espresso_out(filename, index)

    if format == 'nw':
        from ase.io.nwchem import read_nwchem_input
        return read_nwchem_input(filename)

    if format == 'xsd':
        from ase.io.xsd import read_xsd
        return read_xsd(filename)

    raise RuntimeError('File format descriptor ' + format + ' not recognized!')


def write(filename, images, format=None, **kwargs):
    """Write Atoms object(s) to file.

    filename: str
        Name of the file to write to.
    images: Atoms object or list of Atoms objects
        A single Atoms object or a list of Atoms objects.
    format: str
        Used to specify the file-format.  If not given, the
        file-format will be taken from suffix of the filename.

    The accepted output formats:

    =========================  ===========
    format                     short name
    =========================  ===========
    ASE pickle trajectory      traj
    ASE bundle trajectory      bundle
    CUBE file                  cube
    XYZ-file                   xyz
    VASP POSCAR/CONTCAR file   vasp
    ABINIT input file          abinit
    Protein Data Bank          pdb
    CIF-file                   cif
    XCrySDen Structure File    xsf
    FHI-aims geometry file     aims
    gOpenMol .plt file         plt
    Python script              py
    Encapsulated Postscript    eps
    Portable Network Graphics  png
    Persistance of Vision      pov
    VTK XML Image Data         vti
    VTK XML Structured Grid    vts
    VTK XML Unstructured Grid  vtu
    TURBOMOLE coord file       tmol
    exciting                   exi
    AtomEye configuration      cfg
    WIEN2k structure file      struct
    CASTEP cell file           cell
    DftbPlus input file        dftb
    ETSF                       etsf.nc
    DFTBPlus GEN format        gen
    CMR db/cmr-file            db
    CMR db/cmr-file            cmr
    EON reactant.con file      eon
    Gromacs coordinates        gro
    GROMOS96 (only positions)  g96
    X3D                        x3d
    X3DOM HTML                 html
    Extended XYZ file          extxyz
    =========================  ===========

    The use of additional keywords is format specific.

    The ``cube`` and ``plt`` formats accept (plt requires it) a ``data``
    keyword, which can be used to write a 3D array to the file along
    with the nuclei coordinates.

    The ``vti``, ``vts`` and ``vtu`` formats are all specifically directed
    for use with MayaVi, and the latter is designated for visualization of
    the atoms whereas the two others are intended for volume data. Further,
    it should be noted that the ``vti`` format is intended for orthogonal
    unit cells as only the grid-spacing is stored, whereas the ``vts`` format
    additionally stores the coordinates of each grid point, thus making it
    useful for volume date in more general unit cells.

    The ``eps``, ``png``, and ``pov`` formats are all graphics formats,
    and accept the additional keywords:

    rotation: str (default '')
      The rotation angles, e.g. '45x,70y,90z'.

    show_unit_cell: int (default 0)
      Can be 0, 1, 2 to either not show, show, or show all of the unit cell.

    radii: array or float (default 1.0)
      An array of same length as the list of atoms indicating the sphere radii.
      A single float specifies a uniform scaling of the default covalent radii.

    bbox: 4 floats (default None)
      Set the bounding box to (xll, yll, xur, yur) (lower left, upper right).

    colors: array (default None)
      An array of same length as the list of atoms, indicating the rgb color
      code for each atom. Default is the jmol_colors of ase/data/colors.

    scale: int (default 20)
      Number of pixels per Angstrom.

    For the ``pov`` graphics format, ``scale`` should not be specified.
    The elements of the color array can additionally be strings, or 4
    and 5 vectors for named colors, rgb + filter, and rgb + filter + transmit
    specification. This format accepts the additional keywords:

    ``run_povray``, ``display``, ``pause``, ``transparent``,
    ``canvas_width``, ``canvas_height``, ``camera_dist``,
    ``image_plane``, ``camera_type``, ``point_lights``,
    ``area_light``, ``background``, ``textures``, ``celllinewidth``,
    ``bondlinewidth``, ``bondatoms``

    The ``xyz`` format accepts a comment string using the ``comment`` keyword:

    comment: str (default '')
      Optional comment written on the second line of the file.
    """

    if format is None:
        if filename == '-':
            format = 'xyz'
            filename = sys.stdout
        elif 'POSCAR' in filename or 'CONTCAR' in filename:
            format = 'vasp'
        elif 'OUTCAR' in filename:
            format = 'vasp_out'
        elif filename.endswith('etsf.nc'):
            format = 'etsf'
        elif filename.lower().endswith('.con'):
            format = 'eon'
        elif os.path.basename(filename) == 'coord':
            format = 'tmol'
        else:
            suffix = filename.split('.')[-1]
            format = {'cell': 'castep_cell',
                      }.get(suffix, suffix)  # XXX this does not make sense
            # Maybe like this:
##             format = {'traj': 'trajectory',
##                       'nc': 'netcdf',
##                       'exi': 'exciting',
##                       'in': 'aims',
##                       'tmol': 'turbomole',
##                       }.get(suffix, suffix)
            
    if format in ['json', 'db']:
        from ase.db import connect
        connect(filename, format).write(images)
        return
    if format == 'castep_cell':
        from ase.io.castep import write_cell
        write_cell(filename, images, **kwargs)
        return
    if format == 'exi':
        from ase.io.exciting import write_exciting
        write_exciting(filename, images)
        return
    if format == 'cif':
        from ase.io.cif import write_cif
        write_cif(filename, images)
    if format == 'xyz':
        from ase.io.extxyz import write_xyz
        write_xyz(filename, images, columns=['symbols', 'positions'],
                  write_info=False, **kwargs)
        return
    if format == 'extxyz':
        from ase.io.extxyz import write_xyz
        write_xyz(filename, images, **kwargs)
        return
    if format == 'gen':
        from ase.io.gen import write_gen
        write_gen(filename, images)
        return
    elif format == 'in':
        format = 'aims'
    elif format == 'tmol':
        from ase.io.turbomole import write_turbomole
        write_turbomole(filename, images)
        return
    elif format == 'dftb':
        from ase.io.dftb import write_dftb
        write_dftb(filename, images)
        return
    elif format == 'struct':
        from ase.io.wien2k import write_struct
        write_struct(filename, images, **kwargs)
        return
    elif format == 'findsym':
        from ase.io.findsym import write_findsym
        write_findsym(filename, images)
        return
    elif format == 'etsf':
        from ase.io.etsf import ETSFWriter
        writer = ETSFWriter(filename)
        if not isinstance(images, (list, tuple)):
            images = [images]
        writer.write_atoms(images[0])
        writer.close()
        return
    elif format == 'cmr':
        from ase.io.cmr_io import write_db
        return write_db(filename, images, **kwargs)
    elif format == 'eon':
        from ase.io.eon import write_reactant_con
        write_reactant_con(filename, images)
        return
    elif format == 'gro':
        from ase.io.gromacs import write_gromacs
        write_gromacs(filename, images)
        return
    elif format == 'g96':
        from ase.io.gromos import write_gromos
        write_gromos(filename, images)
        return
    elif format == 'html':
        from ase.io.x3d import write_html
        write_html(filename, images)
        return

    format = {'traj': 'trajectory',
              'nc': 'netcdf',
              'bundle': 'bundletrajectory'
              }.get(format, format)
    name = 'write_' + format

    if format in ['vti', 'vts', 'vtu']:
        format = 'vtkxml'

    if format is None:
        format = filetype(filename)

    try:
        write = getattr(__import__('ase.io.%s' % format, {}, {}, [name]), name)
    except ImportError:
        raise TypeError('Unknown format: "%s".' % format)

    write(filename, images, **kwargs)


def string2index(string):
    if ':' not in string:
        return int(string)
    i = []
    for s in string.split(':'):
        if s == '':
            i.append(None)
        else:
            i.append(int(s))
    i += (3 - len(i)) * [None]
    return slice(*i)


def filetype(filename):
    """Try to guess the type of the file."""
    if os.path.isdir(filename):
        # Potentially a BundleTrajectory
        if BundleTrajectory.is_bundle(filename):
            return 'bundle'
        elif os.path.normpath(filename) == 'states':
            return 'eon'
        else:
            raise IOError('Directory: ' + filename)

    if filename.startswith('pg://'):
        return 'postgresql'

    fileobj = open(filename, 'rU')
    s3 = fileobj.read(3)
    if len(s3) == 0:
        raise IOError('Empty file: ' + filename)

    if s3.startswith('{"'):
        return 'json'

    if filename.endswith('.db'):
        return 'db'

    if filename.lower().endswith('.cmr'):
        return 'cmr'

    if is_tarfile(filename):
        return 'gpw'

    if s3 == 'CDF':
        from ase.io.pupynere import NetCDFFile
        nc = NetCDFFile(filename)
        if 'number_of_dynamic_atoms' in nc.dimensions:
            return 'dacapo'

        history = nc.history
        if history == 'GPAW restart file':
            return 'gpw-nc'
        if history == 'ASE trajectory':
            return 'nc'
        if history == 'Dacapo':
            return 'dacapo'
        if hasattr(nc, 'file_format') and nc.file_format.startswith('ETSF'):
            return 'etsf'
        raise IOError('Unknown netCDF file!')

    if is_zipfile(filename):
        return 'vnl'

    fileobj.seek(0)
    lines = fileobj.readlines(1000)

    if lines[0].startswith('PickleTrajectory'):
        return 'traj'

    if (lines[1].startswith('OUTER LOOP:') or
        filename.lower().endswith('.cube')):
        return 'cube'

    if '  ___ ___ ___ _ _ _  \n' in lines:
        return 'gpaw-text'

    if (' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n'
        in lines[:90]):
        return 'dacapo-text'

    for line in lines:
        if line[0] != '#':
            word = line.strip()
            if word in ['ANIMSTEPS', 'CRYSTAL', 'SLAB', 'POLYMER', 'MOLECULE']:
                return 'xsf'

    filename_v = os.path.basename(filename)
    if 'POSCAR' in filename_v or 'CONTCAR' in filename_v:
        return 'vasp'

    if 'OUTCAR' in filename_v:
        return 'vasp_out'

    if 'XDATCAR' in filename_v:
        return 'vasp_xdatcar'

    if filename.lower().endswith('.exi'):
        return 'exi'

    if filename.lower().endswith('.mol'):
        return 'mol'

    if filename.lower().endswith('.pdb'):
        return 'pdb'

    if filename.lower().endswith('.cif'):
        return 'cif'

    if filename.lower().endswith('.struct'):
        return 'struct'

    if filename.lower().endswith('.struct_out'):
        return 'struct_out'

    fileobj.seek(0)
    while True:
        line = fileobj.readline()
        if not line:
            break
        if 'Invoking FHI-aims ...' in line:
            return 'aims_out'
        if 'atom' in line:
            data = line.split()
            try:
                Atoms(symbols=[data[4]],
                      positions=[[float(data[1]),
                                  float(data[2]),
                                  float(data[3])]])
                return 'aims'
            except:
                pass

    if filename.lower().endswith('.in'):
        fileobj.seek(0)
        while True:
            line = fileobj.readline()
            if not line:
                break
            if ('&system' in line) or ('&SYSTEM' in line):
                return 'esp_in'
        return 'aims'

    if filename.lower().endswith('.cfg'):
        return 'cfg'

    if os.path.split(filename)[1] == 'atoms.dat':
        return 'iwm'

    if filename.endswith('I_info'):
        return 'Cmdft'

    if lines[0].startswith('$coord') or os.path.basename(filename) == 'coord':
        return 'tmol'

    if (lines[0].startswith('$grad') or
        os.path.basename(filename) == 'gradient'):
        return 'tmol-gradient'

    if lines[0].startswith('Geometry'):
        return 'dftb'

    if filename.lower().endswith('.geom'):
        return 'castep_geom'

    if filename.lower().endswith('.castep'):
        return 'castep'

    if filename.lower().endswith('.cell'):
        return 'castep_cell'
    if s3 == '<?x' and not filename.endswith('xsd'):
        from ase.io.vtkxml import probe_vtkxml
        xmltype = probe_vtkxml(filename)
        if xmltype == 'ImageData':
            return 'vti'
        elif xmltype == 'StructuredGrid':
            return 'vts'
        elif xmltype == 'UnstructuredGrid':
            return 'vtu'
        elif xmltype is not None:
            raise IOError('Unknown VTK XML file!')

    if filename.lower().endswith('.sdf'):
        return 'sdf'

    if filename.lower().endswith('.gen'):
        return 'gen'

    if filename.lower().endswith('.con'):
        return 'eon'

    if 'ITEM: TIMESTEP\n' in lines:
        return 'lammps'

    if filename.lower().endswith('.gro'):
        return 'gromacs'

    if filename.lower().endswith('.log'):
        return 'gaussian_out'

    if filename.lower().endswith('.com'):
        return 'gaussian'

    if filename.lower().endswith('.g96'):
        return 'gromos'

    if filename.lower().endswith('.out'):
        return 'esp_out'

    if filename.endswith('.nw'):
        return 'nw'

    if filename.endswith('xsd'):
        return 'xsd'

    return 'xyz'
