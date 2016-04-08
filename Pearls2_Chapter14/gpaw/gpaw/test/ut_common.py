"""Common code base for maintaining backwards compatibility in ut_xxx tests."""

__all__ = ['ase_svnversion', 'shapeopt', 'TestCase', 'TextTestRunner', \
    'CustomTextTestRunner', 'defaultTestLoader', 'initialTestLoader', \
    'create_random_atoms', 'create_parsize_maxbands', 'create_parsize_minbands']

partest = True

# -------------------------------------------------------------------

# Maintain backwards compatibility with ASE 3.1.0 svn. rev. 1158 or later
try:
    from ase.svnversion import svnversion as ase_svnversion
except ImportError:
    # Fall back on minimum required ASE svn.rev.
    ase_svnversion = 1158
else:
    # From test/ase3k_version.py.
    full_ase_svnversion = ase_svnversion
    if ase_svnversion[-1] == 'M':
        ase_svnversion = ase_svnversion[:-1]
    if ase_svnversion.rfind(':') != -1:
        ase_svnversion = ase_svnversion[:ase_svnversion.rfind(':')]
    ase_svnversion = int(ase_svnversion)

# Using a feature from ASE 3.1.0 svn. rev. 1001 or later.
from ase.utils.memory import shapeopt

if partest:
    from gpaw.test.parunittest import ParallelTestCase as TestCase, \
        ParallelTextTestRunner as TextTestRunner, ParallelTextTestRunner as \
        CustomTextTestRunner, defaultParallelTestLoader as defaultTestLoader
    def CustomTextTestRunner(logname, verbosity=1):
        return TextTestRunner(stream=logname, verbosity=verbosity)
else:
    # Using a features from ASE 3.1.0 svn. rev. 929 or later.
    from ase.test import CustomTestCase as TestCase, CustomTextTestRunner
    from unittest import TextTestRunner, defaultTestLoader

from copy import copy
initialTestLoader = copy(defaultTestLoader)
assert hasattr(initialTestLoader, 'testMethodPrefix')
initialTestLoader.testMethodPrefix = 'verify'

# -------------------------------------------------------------------

import numpy as np

from math import sin, cos
from ase import Atoms
from ase.structure import molecule
from ase.units import Bohr
from gpaw.mpi import compare_atoms
from gpaw.utilities.tools import md5_array

def create_random_atoms(gd, nmolecules=10, name='NH2', mindist=4.5 / Bohr):
    """Create gas-like collection of atoms from randomly placed molecules.
    Applies rigid motions to molecules, translating the COM and/or rotating
    by a given angle around an axis of rotation through the new COM. These
    atomic positions obey the minimum distance requirement to zero-boundaries.

    Warning: This is only intended for testing parallel grid/LFC consistency.
    """
    atoms = Atoms(cell=gd.cell_cv * Bohr, pbc=gd.pbc_c)

    # Store the original state of the random number generator
    randstate = np.random.get_state()
    seed = np.array([md5_array(data, numeric=True)
                     for data in
                     [nmolecules, gd.cell_cv, gd.pbc_c, gd.N_c]]).astype(int)
    np.random.seed(seed % 4294967296)

    for m in range(nmolecules):
        amol = molecule(name)
        amol.set_cell(gd.cell_cv * Bohr)

        # Rotate the molecule around COM according to three random angles
        # The rotation axis is given by spherical angles phi and theta
        v,phi,theta = np.random.uniform(0.0, 2*np.pi, 3) # theta [0,pi[ really
        axis = np.array([cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)])
        amol.rotate(axis, v)

        # Find the scaled length we must transverse along the given axes such
        # that the resulting displacement vector is `mindist` from the cell
        # face corresponding to that direction (plane with unit normal n_v).
        sdist_c = np.empty(3)
        if not gd.orthogonal:
            for c in range(3):
                n_v = gd.xxxiucell_cv[c] / np.linalg.norm(gd.xxxiucell_cv[c])
                sdist_c[c] = mindist / np.dot(gd.cell_cv[c], n_v)
        else:
            sdist_c[:] = mindist / gd.cell_cv.diagonal()
        assert np.all(sdist_c > 0), 'Displacment vectors must be inside cell.'

        # Scaled dimensions of the smallest possible box centered on the COM
        spos_ac = amol.get_scaled_positions() # NB! must not do a "% 1.0"
        scom_c = np.dot(gd.icell_cv, amol.get_center_of_mass())
        sbox_c = np.abs(spos_ac-scom_c[np.newaxis,:]).max(axis=0)
        sdelta_c = (1-np.array(gd.pbc_c)) * (sbox_c + sdist_c)
        assert (sdelta_c < 1.0-sdelta_c).all(), 'Box is too tight to fit atoms.'
        scenter_c = [np.random.uniform(d,1-d) for d in sdelta_c]
        center_v = np.dot(scenter_c, gd.cell_cv)

        # Translate the molecule such that COM is located at random center
        offset_av = (center_v-amol.get_center_of_mass()/Bohr)[np.newaxis,:]
        amol.set_positions(amol.get_positions()+offset_av*Bohr)
        assert np.linalg.norm(center_v-amol.get_center_of_mass()/Bohr) < 1e-9
        atoms.extend(amol)

    # Restore the original state of the random number generator
    np.random.set_state(randstate)
    assert compare_atoms(atoms)
    return atoms


# -------------------------------------------------------------------

from gpaw.utilities import gcd
from gpaw import parsize_domain, parsize_bands

def create_parsize_maxbands(nbands, world_size):
    """Safely parse command line parallel arguments for band parallel case."""
    # D: number of domains
    # B: number of band groups   
    if parsize_bands is None:
        if parsize_domain is None:
            B = gcd(nbands, world_size) # largest possible
            D = world_size // B
        else:
            D = parsize_domain
            B = gcd(nbands, world_size // np.prod(D))
    else:
        B = parsize_bands
        D = parsize_domain or world_size // B
    return D, B

def create_parsize_minbands(nbands, world_size):
    __doc__ = create_parsize_maxbands.__doc__
    return create_parsize_maxbands(1, world_size)
