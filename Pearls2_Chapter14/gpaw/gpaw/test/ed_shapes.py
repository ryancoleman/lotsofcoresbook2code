from ase import Atoms
from gpaw.fdtd.poisson_fdtd import QSFDTD
from gpaw.fdtd.polarizable_material import PermittivityPlus, PolarizableMaterial, \
                                           PolarizableSphere, PolarizableBox, \
                                           PolarizableEllipsoid, PolarizableRod, \
                                           PolarizableTetrahedron
from gpaw.mpi import world
from gpaw.tddft import photoabsorption_spectrum, units
from gpaw.test import equal
import numpy as np

# Whole simulation cell (Angstroms)
cell = [40, 40, 20];

# Classical subsystem
classical_material = PolarizableMaterial()
classical_material.add_component(PolarizableSphere(permittivity = PermittivityPlus(data = [[1.20, 0.20, 25.0]]),
                                                   center = [10, 10, 10],
                                                   radius = 4.5))
classical_material.add_component(PolarizableBox(permittivity = PermittivityPlus(data = [[1.40, 0.20, 25.0]]),
                                                corner1 = [18.1, 5.1, 5.1],
                                                corner2 = [22.9, 14.9, 14.9]))
classical_material.add_component(PolarizableEllipsoid(permittivity = PermittivityPlus(data = [[1.60, 0.20, 25.0]]),
                                                center  = [30.0, 10.0, 10.0],
                                                radii   = [ 3.9, 5.9, 4.9]))
classical_material.add_component(PolarizableRod(permittivity = PermittivityPlus(data = [[1.80, 0.20, 25.0]]),
                                                corners = [[10.0, 21.5, 10.0], [10.0, 33.5, 10.0]],
                                                round_corners = True,
                                                radius = 3.9))
classical_material.add_component(PolarizableRod(permittivity = PermittivityPlus(data = [[1.00, 0.20, 25.0]]),
                                                corners = [[20.0, 21.5, 10.0], [25.0, 33.5, 10.0]],
                                                round_corners = False,
                                                radius = 2.9))
classical_material.add_component(PolarizableTetrahedron(permittivity = PermittivityPlus(data = [[0.80, 0.20, 25.0]]),
                                                        corners = [[24.1, 16.1, 5.1],
                                                                   [30.1, 36.1, 6.1],
                                                                   [36.4, 27.6, 7.1],
                                                                   [30.0, 25.0, 14.9]]))


# Wrap calculators
qsfdtd = QSFDTD(classical_material = classical_material,
                atoms              = None,
                cells              = (cell, 2.00),
                spacings           = [1.60, 0.40],
                remove_moments     = (1, 1))

# Run
energy = qsfdtd.ground_state('gs.gpw', eigensolver = 'cg', nbands = -1)
qsfdtd.time_propagation('gs.gpw', kick_strength=[0.000, 0.000, 0.001], time_step=10, iterations=5, dipole_moment_file='dmCl.dat')

# Restart and run
qsfdtd.write('td.gpw', mode='all')
qsfdtd.time_propagation('td.gpw', kick_strength=None, time_step=10, iterations=5, dipole_moment_file='dmCl.dat')

# Test
ref_cl_dipole_moment = [-1.01218372e-04, -3.03603352e-05, 1.86063694e-01]
tol = 0.0001
equal(qsfdtd.td_calc.hamiltonian.poisson.get_classical_dipole_moment(), ref_cl_dipole_moment, tol)

