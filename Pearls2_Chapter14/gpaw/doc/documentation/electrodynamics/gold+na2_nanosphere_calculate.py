from ase import Atoms
from ase.units import Hartree
from gpaw.fdtd.poisson_fdtd import QSFDTD
from gpaw.fdtd.polarizable_material import PermittivityPlus, PolarizableMaterial, PolarizableSphere
from gpaw.tddft import photoabsorption_spectrum
import numpy as np

# Nanosphere radius (Angstroms)
radius = 7.40

# Geometry
atom_center     = np.array([30., 15., 15.])
sphere_center   = np.array([15., 15., 15.])
simulation_cell = np.array([40., 30., 30.])

# Atoms object
atoms = Atoms('Na2', atom_center + np.array([[-1.5, 0.0, 0.0],
                                             [ 1.5, 0.0, 0.0]]))

# Permittivity of Gold (from  J. Chem. Phys. 137, 074113 (2012))
eps_gold = PermittivityPlus(data = [[0.2350, 0.1551,  95.62],
                                    [0.4411, 0.1480, -12.55],
                                    [0.7603,  1.946, -40.89],
                                    [1.161,   1.396,  17.22],
                                    [2.946,   1.183,  15.76],
                                    [4.161,   1.964,  36.63],
                                    [5.747,   1.958,  22.55],
                                    [7.912,   1.361,  81.04]])

# 1) Nanosphere only
classical_material = PolarizableMaterial()                            

classical_material.add_component(PolarizableSphere(center = sphere_center,
                                 radius                   = radius,
                                 permittivity             = eps_gold))

qsfdtd = QSFDTD(classical_material = classical_material,
                atoms              = None,
                cells              = simulation_cell,
                spacings           = [2.0, 0.5],
                remove_moments     = (1, 1))

energy = qsfdtd.ground_state('gs.gpw',
                             nbands = 1)

qsfdtd.time_propagation('gs.gpw',
                        kick_strength=[0.001, 0.000, 0.000],
                        time_step=10,
                        iterations=1500,
                        dipole_moment_file='dm.dat')

photoabsorption_spectrum('dm.dat', 'spec.1.dat', width=0.15)


# 2) Na2 only (radius=0)
classical_material = PolarizableMaterial()                            
classical_material.add_component(PolarizableSphere(center = sphere_center,
                                 radius                   = 0.0,
                                 permittivity             = eps_gold))

qsfdtd = QSFDTD(classical_material = classical_material,
                atoms              = atoms,
                cells              = (simulation_cell, 4.0), # vacuum = 4.0 Ang
                spacings           = [2.0, 0.5],
                remove_moments     = (1, 1))

energy = qsfdtd.ground_state('gs.gpw', nbands = -1)

qsfdtd.time_propagation('gs.gpw',
                        kick_strength=[0.001, 0.000, 0.000],
                        time_step=10,
                        iterations=1500,
                        dipole_moment_file='dm.dat')

photoabsorption_spectrum('dm.dat', 'spec.2.dat', width=0.15)

# 3) Nanosphere + Na2
classical_material = PolarizableMaterial()                            
classical_material.add_component(PolarizableSphere(center = sphere_center,
                                 radius                   = radius,
                                 permittivity             = eps_gold))

qsfdtd = QSFDTD(classical_material = classical_material,
                atoms              = atoms,
                cells              = (simulation_cell, 4.0), # vacuum = 4.0 Ang
                spacings           = [2.0, 0.5],
                remove_moments     = (1, 1))

energy = qsfdtd.ground_state('gs.gpw', nbands = -1)

qsfdtd.time_propagation('gs.gpw',
                        kick_strength=[0.001, 0.000, 0.000],
                        time_step=10,
                        iterations=1500,
                        dipole_moment_file='dm.dat')

photoabsorption_spectrum('dm.dat', 'spec.3.dat', width=0.15)
