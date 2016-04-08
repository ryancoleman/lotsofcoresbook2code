from ase import Atoms
from gpaw import GPAW
from gpaw.fdtd.poisson_fdtd import FDTDPoissonSolver
from gpaw.fdtd.polarizable_material import PermittivityPlus, PolarizableMaterial, PolarizableSphere
from gpaw.mpi import world
from gpaw.tddft import TDDFT, photoabsorption_spectrum, units
from gpaw.test import equal
import numpy as np

# Whole simulation cell (Angstroms)
large_cell = [20, 20, 30];

# Quantum subsystem
atom_center = np.array([10.0, 10.0, 20.0]);
atoms = Atoms('Na2', [atom_center + [0.0, 0.0, -1.50],
                      atom_center + [0.0, 0.0, +1.50]]);

# Permittivity file
if world.rank==0:
    fo = open("ed.txt", "wb")
    fo.writelines(["1.20 0.20 25.0"])
    fo.close()
world.barrier()

# Classical subsystem
classical_material = PolarizableMaterial()
sphere_center = np.array([10.0, 10.0, 10.0]);
classical_material.add_component(PolarizableSphere(permittivity = PermittivityPlus('ed.txt'),
                                                   center = sphere_center,
                                                   radius = 5.0
                                                   ))

# Combined Poisson solver
poissonsolver = FDTDPoissonSolver(classical_material  = classical_material,
                                  qm_spacing          = 0.40,
                                  cl_spacing          = 0.40*4,
                                  cell                = large_cell,
                                  remove_moments      = (1, 4),
                                  communicator        = world,
                                  potential_coupler   = 'Refiner')
poissonsolver.set_calculation_mode('iterate')

# Combined system
atoms.set_cell(large_cell)
atoms, qm_spacing, gpts = poissonsolver.cut_cell(atoms,
                                                 vacuum=2.50)

# Initialize GPAW
gs_calc = GPAW(gpts          = gpts,
               eigensolver   = 'cg',
               nbands        = -1,
               poissonsolver = poissonsolver);
atoms.set_calculator(gs_calc)

# Ground state
energy = atoms.get_potential_energy()

# Save state
gs_calc.write('gs.gpw', 'all')
classical_material = None
gs_calc = None

# Initialize TDDFT and FDTD
kick = [0.0, 0.0, 1.0e-3]
time_step = 10.0
max_time = 100 # 0.1 fs

td_calc = TDDFT('gs.gpw')
td_calc.absorption_kick(kick_strength=kick)
td_calc.hamiltonian.poisson.set_kick(kick)

# Propagate TDDFT and FDTD
td_calc.propagate(time_step,  max_time/time_step/2, 'dm.dat', 'td.gpw')

td_calc2 = TDDFT('td.gpw')
td_calc2.propagate(time_step,  max_time/time_step/2, 'dm.dat', 'td.gpw')

# Test
ref_cl_dipole_moment = [ -5.16149623e-14, -5.89090408e-14,  3.08450150e-02]
ref_qm_dipole_moment = [ -2.63340461e-11,  2.61812794e-12, -9.35619772e-02]
tol = 0.0001
equal(td_calc2.hamiltonian.poisson.get_classical_dipole_moment(), ref_cl_dipole_moment, tol)
equal(td_calc2.hamiltonian.poisson.get_quantum_dipole_moment(), ref_qm_dipole_moment, tol)

