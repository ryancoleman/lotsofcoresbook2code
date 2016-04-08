from ase import Atoms
from ase.units import Hartree
from gpaw.fdtd.poisson_fdtd import QSFDTD
from gpaw.fdtd.polarizable_material import PermittivityPlus, PolarizableMaterial, PolarizableSphere
from gpaw.tddft import photoabsorption_spectrum
from gpaw.mpi import world, serial_comm
import numpy as np

# Nanosphere radius (Angstroms)
radius = 50.0

# Whole simulation cell (Angstroms)
large_cell = np.array([3.*radius, 3.*radius, 3.*radius])

# Permittivity of Gold (from  J. Chem. Phys. 137, 074113 (2012))
gold = [[0.2350, 0.1551,  95.62],
        [0.4411, 0.1480, -12.55],
        [0.7603,  1.946, -40.89],
        [1.161,   1.396,  17.22],
        [2.946,   1.183,  15.76],
        [4.161,   1.964,  36.63],
        [5.747,   1.958,  22.55],
        [7.912,   1.361,  81.04]]

# Initialize classical material
classical_material = PolarizableMaterial()

# Classical nanosphere
classical_material.add_component(
        PolarizableSphere(center = 0.5*large_cell,
                          radius = radius,
                          permittivity = PermittivityPlus(data=gold))
        )

# Quasistatic FDTD
qsfdtd = QSFDTD(classical_material = classical_material,
                atoms              = None,
                cells              = large_cell,
                spacings           = [8.0, 1.0],
                communicator       = world,
                remove_moments     = (4, 1))
# Run ground state
energy = qsfdtd.ground_state('gs.gpw', nbands = -1)

# Run time evolution
qsfdtd.time_propagation('gs.gpw',
                        time_step=10,
                        iterations=1000,
                        kick_strength=[0.001, 0.000, 0.000],
                        dipole_moment_file='dm.dat')
# Spectrum
photoabsorption_spectrum('dm.dat', 'spec.dat', width=0.0)


