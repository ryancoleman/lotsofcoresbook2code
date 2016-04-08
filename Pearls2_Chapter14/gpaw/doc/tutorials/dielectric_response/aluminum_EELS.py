import numpy as np
from ase.lattice import bulk
from gpaw import GPAW
from gpaw.response.df import DielectricFunction

# Part 1: Ground state calculation
atoms = bulk('Al', 'fcc', a=4.043)      # Generate fcc crystal structure for aluminum
calc = GPAW(mode='pw',                  # GPAW calculator initialization
            kpts={'density': 5.0, 'gamma': True})

atoms.set_calculator(calc)
atoms.get_potential_energy()            # Ground state calculation is performed
calc.write('Al.gpw','all')              # Use 'all' option to write wavefunctions

# Part 2: Spectrum calculation
df = DielectricFunction(calc='Al.gpw')  # Ground state gpw file as input

q_c = [1.0 / 13, 0, 0]          # Momentum transfer, must be the difference between two kpoints !
df.get_eels_spectrum(q_c=q_c)           # By default, a file called 'eels.csv' is generated
