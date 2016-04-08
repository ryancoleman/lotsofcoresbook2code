#!/usr/bin/env python

import os

from ase import Atoms
from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton
from ase.io import write
from ase.structure import molecule
from ase.md.verlet import VelocityVerlet
from ase.md import MDLogger
from ase.units import *
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.dftb import read_dftb_velocities, write_dftb_velocities
from ase.io import read, write

o2 = molecule('O2')
h2_1 = molecule('H2')
h2_2 = molecule('H2')
o2.translate([0,0.01,0])
h2_1.translate([0,0,3])
h2_1.rotate_euler(center='COP',theta=3.1415/2)
h2_2.translate([0,0,-3])
h2_2.rotate_euler(center='COP',theta=3.1415/2)
o2.set_velocities(([0,0,0],[0,0,0]))
h2_1.set_velocities(([0,0,-3.00],[0,0,-3.000]))
h2_2.set_velocities(([0,0,3.000],[0,0,3.000]))
test = o2 + h2_1 + h2_2

#1fs = 41.3 au
#1000K = 0.0031668 au
calculator_NVE = Dftb(label='h2o',atoms=test,
run_manyDftb_steps = True,
Hamiltonian_MaxAngularMomentum_ = '',
Hamiltonian_MaxAngularMomentum_O = '"p"',
Hamiltonian_MaxAngularMomentum_H = '"s"',
Driver_ = 'VelocityVerlet',
Driver_MDRestartFrequency = 10,
Driver_Velocities_ = '',
Driver_Velocities_empty = '<<+ "velocities.txt"',
Driver_Steps = 1000, 
Driver_KeepStationary = 'Yes', 
Driver_TimeStep = 4.13, 
Driver_Thermostat_ = 'None',
Driver_Thermostat_empty = '',
)

#1fs = 41.3 au
#1000K = 0.0031668 au
calculator_NVT = Dftb(label='h2o',atoms=test,
run_manyDftb_steps = True,
Hamiltonian_MaxAngularMomentum_ = '',
Hamiltonian_MaxAngularMomentum_O = '"p"',
Hamiltonian_MaxAngularMomentum_H = '"s"',
Driver_ = 'VelocityVerlet',
Driver_MDRestartFrequency = 5,
Driver_Velocities_ = '',
Driver_Velocities_empty = '<<+ "velocities.txt"',
Driver_Steps = 500, 
Driver_KeepStationary = 'Yes', 
Driver_TimeStep = 8.26, 
Driver_Thermostat_ = 'Berendsen',
#Driver_Thermostat_Temperature = 0.00339845142, # 800 deg Celcius
Driver_Thermostat_Temperature = 0.0, # 0 deg Kelvin
Driver_Thermostat_CouplingStrength = 0.01,
)

write_dftb_velocities(test, 'velocities.txt')
os.system("rm md.log.* md.out* geo_end*xyz")
test.set_calculator(calculator_NVE)
dyn = VelocityVerlet(test, 0.000*fs)  #  fs time step.
dyn.attach(MDLogger(dyn, test, 'md.log.NVE', header=True, stress=False,
                    peratom=False, mode="w"), interval=1)
dyn.run(1) # run NVE ensemble using DFTB's own driver
test = read('geo_end.gen')
write('test.afterNVE.xyz', test)

read_dftb_velocities(test, filename='geo_end.xyz')
write_dftb_velocities(test, 'velocities.txt')

os.system("mv md.out md.out.NVE")
os.system("mv geo_end.xyz geo_end_NVE.xyz")

test.set_calculator(calculator_NVT)
os.system("rm md.log.NVT")
dyn.attach(MDLogger(dyn, test, 'md.log.NVT', header=True, stress=False,
                    peratom=False, mode="w"), interval=1)
dyn.run(1)  # run NVT ensemble using DFTB's own driver
test = read('geo_end.gen')
read_dftb_velocities(test, filename='geo_end.xyz')

os.system("mv md.out md.out.NVT")
os.system("mv geo_end.xyz geo_end_NVT.xyz")
write_dftb_velocities(test, 'velocities.txt')

# to get a movie of the full collision:
# cat geo_end_NVE.xyz geo_end_NVT.xyz > all.xyz
# and:
#  ag all.xyz
