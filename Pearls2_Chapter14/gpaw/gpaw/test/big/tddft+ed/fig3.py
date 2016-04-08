from ase import Atoms
from ase.io import read
from gpaw import GPAW, Mixer
from gpaw.fdtd import FDTDPoissonSolver, PermittivityPlus, PolarizableMaterial, PolarizableSphere
from gpaw.mpi import world
from gpaw.svnversion import svnversion
from gpaw.tddft import TDDFT, photoabsorption_spectrum, units
from gpaw.test import equal
import numpy as np
import sys

# Command line arguments
kwargs         = dict(x.split('=', 1) for x in sys.argv[1:])
coupler        = kwargs.pop('coupler', 'Refiner')
coupling_level = kwargs.pop('level', 'both')
maxL_qm        = int(kwargs.pop('maxL_qm', 4))
maxL_cl        = int(kwargs.pop('maxL_cl', 1))
radius         = float(kwargs.pop('radius'  , 5.0))
distance       = float(kwargs.pop('distance', 10.0))
qm_spacing     = float(kwargs.pop('qmh', 0.50))
cl_spacing     = float(kwargs.pop('clh', 2.00))
qm_vacuum      = float(kwargs.pop('qm_vacuum', 5.00))
cx             = float(kwargs.pop('cx', 96.0))
cy             = float(kwargs.pop('cy', 48.0))
cz             = float(kwargs.pop('cz', 48.0))
time_step      = float(kwargs.pop('dt', 20))
title          = kwargs.pop('title', 'NP2+Na2')
kick           = np.array([float(kwargs.pop('kick', 0.001)), 0.0, 0.0])
npworld = world.size
bondlength    = 3.0788

# Permittivity
if world.rank==0:
    fo = open("silver.txt", "wb")
    fo.writelines(["0.1696  0.1795 135.0\n",
                   "0.3655  0.2502 -40.30\n",
                   "0.6312  2.114  -50.06\n",
                   "1.175   1.627   16.73\n",
                   "2.077   1.820    7.651\n",
                   "4.018   1.049  -15.36\n",
                   "4.243   0.9967  18.07\n",
                   "5.303   2.592   40.42\n",
                   "7.197   2.774   31.02"]);
    fo.close()

# Whole simulation cell (Angstroms)
large_cell = [cx, cy, cz];

# Quantum subsystem
atom_center = np.array([0.5*cx, 0.5*cy, 0.5*cz]);
atoms = Atoms('Na2', [[-0.5*bondlength, 0.0, 0.0],
                      [ 0.5*bondlength, 0.0, 0.0]])

atoms.center(vacuum=0.0)
atoms.positions = atom_center + (atoms.positions - 0.5*np.diag(atoms.get_cell()))

classical_material = PolarizableMaterial()
tag = 'dummy'

# Classical subsystem
classical_material.add_component(PolarizableSphere(center = [0.5*cx-0.5*distance-radius, 0.5*cy, 0.5*cz],
                                                   radius = radius,
                                                   permittivity = PermittivityPlus('silver.txt')))
classical_material.add_component(PolarizableSphere(center = [0.5*cx+0.5*distance+radius, 0.5*cy, 0.5*cz],
                                                   radius = radius,
                                                   permittivity = PermittivityPlus('silver.txt')))

# Combined Poisson solver
poissonsolver = FDTDPoissonSolver(classical_material  = classical_material,
                                  qm_spacing          = qm_spacing,
                                  cl_spacing          = cl_spacing,
                                  cell                = large_cell,
                                  remove_moments      = (maxL_qm, maxL_cl),
                                  communicator        = world,
                                  debug_plots         = 0,
                                  potential_coupler   = coupler,
                                  coupling_level      = coupling_level,
                                  dm_fname            = 'dmCl.%s.dat' % tag,
                                  tag                 = tag)
poissonsolver.set_calculation_mode('iterate')

# Combined system
atoms.set_cell(large_cell)
atoms, qm_spacing, gpts = poissonsolver.cut_cell(atoms,
                                                 vacuum=qm_vacuum)

# Initialize GPAW
gs_calc = GPAW(gpts          = gpts,
               eigensolver   = 'cg',
               nbands        = -2,
               poissonsolver = poissonsolver)
atoms.set_calculator(gs_calc)

# Ground state
energy = atoms.get_potential_energy()

# Save state
gs_calc.write('gs.%s.gpw' % tag, 'all')

# Initialize TDDFT and FDTD
td_calc = TDDFT('gs.%s.gpw' % tag)
td_calc.absorption_kick(kick_strength=kick)
td_calc.hamiltonian.poisson.set_kick(kick)

# Propagate TDDFT and FDTD
td_calc.propagate(time_step,
                  50,
                  'dm.%s.dat' % tag,
                  'td.%s.gpw' % tag)

# Test
ref_cl_dipole_moment = [ 2.72623607e-02,  1.98393701e-09, -1.98271199e-09 ]
ref_qm_dipole_moment = [ 1.44266213e-02,  1.04985435e-09, -1.04920610e-09 ]

tol = 0.0001
equal(td_calc.get_dipole_moment(), ref_qm_dipole_moment, tol)
equal(td_calc.hamiltonian.poisson.get_dipole_moment(), ref_cl_dipole_moment, tol)

