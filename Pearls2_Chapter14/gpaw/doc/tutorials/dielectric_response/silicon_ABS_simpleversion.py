import numpy as np
from ase.lattice import bulk
from gpaw import GPAW
from gpaw.response.df import DielectricFunction

# Part 1: Ground state calculation
atoms = bulk('Si', 'diamond', a=5.431)   # Generate diamond crystal structure for silicon
calc = GPAW(mode='pw', kpts=(4,4,4))     # GPAW calculator initialization
 
atoms.set_calculator(calc)               
atoms.get_potential_energy()             # Ground state calculation is performed
calc.write('si.gpw', 'all')              # Use 'all' option to write wavefunction

# Part 2 : Spectrum calculation          # DF: dielectric function object
df = DielectricFunction(calc='si.gpw',   # Ground state gpw file (with wavefunction) as input
                        domega0=0.05)    # Using nonlinear frequency grid
df.get_dielectric_function()             # By default, a file called 'df.csv' is generated

