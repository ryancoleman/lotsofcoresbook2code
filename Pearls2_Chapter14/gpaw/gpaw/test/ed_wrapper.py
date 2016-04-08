from ase import Atoms
from gpaw.fdtd.poisson_fdtd import QSFDTD
from gpaw.fdtd.polarizable_material import PermittivityPlus, PolarizableMaterial, PolarizableSphere
from gpaw.mpi import world
from gpaw.tddft import photoabsorption_spectrum, units
from gpaw.test import equal
import numpy as np

# Whole simulation cell (Angstroms)
cell = [20, 20, 30];

# Quantum subsystem
atom_center = np.array([10.0, 10.0, 20.0]);
atoms = Atoms('Na2', [atom_center + [0.0, 0.0, -1.50],
                      atom_center + [0.0, 0.0, +1.50]]);

# Classical subsystem
sphere_center = np.array([10.0, 10.0, 10.0]);
classical_material = PolarizableMaterial()
classical_material.add_component(PolarizableSphere(permittivity = PermittivityPlus(data = [[1.20, 0.20, 25.0]]),
                                                   center = sphere_center,
                                                   radius = 5.0
                                                   ))

# Wrap calculators
qsfdtd = QSFDTD(classical_material = classical_material,
                atoms              = atoms,
                cells              = (cell, 2.50),
                spacings           = [1.60, 0.40],
                remove_moments     = (1, 4),
                communicator       = world)

# Run
energy = qsfdtd.ground_state('gs.gpw', eigensolver = 'cg', nbands = -1)
qsfdtd.time_propagation('gs.gpw', kick_strength=[0.000, 0.000, 0.001], time_step=10, iterations=5, dipole_moment_file='dm.dat', restart_file='td.gpw')
qsfdtd.time_propagation('td.gpw', kick_strength=None, time_step=10, iterations=5, dipole_moment_file='dm.dat')

# Test
ref_cl_dipole_moment = [ -5.16149623e-14, -5.89090408e-14,  3.08450150e-02]
ref_qm_dipole_moment = [ -2.63340461e-11,  2.61812794e-12, -9.35619772e-02]
tol = 0.0001
equal(qsfdtd.td_calc.hamiltonian.poisson.get_classical_dipole_moment(), ref_cl_dipole_moment, tol)
equal(qsfdtd.td_calc.hamiltonian.poisson.get_quantum_dipole_moment(), ref_qm_dipole_moment, tol)

