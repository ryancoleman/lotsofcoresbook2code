import os
from subprocess import Popen, PIPE, STDOUT

from ase import Atoms
from ase.calculators.turbomole import Turbomole

# Delete old coord, control, ... files, if exist
for f in ['coord',
          'basis',
          'energy',
          'mos',
          'statistics',
          'control']:
    if os.path.exists(f):
        os.remove(f)

atoms = Atoms('H2', positions=[(0, 0, 0), (0, 0, 1.1)])
atoms.set_calculator(Turbomole()) # Writes a coord file as well

# Write all commands for the define command in a string
define_str = '\n\na coord\n*\nno\nb all sto-3g hondo\n*\neht\n\n\n\n*'

# Run define
p = Popen('define', stdout=PIPE, stdin=PIPE, stderr=STDOUT)
stdout = p.communicate(input=define_str)

# Run turbomole
atoms.get_potential_energy()
