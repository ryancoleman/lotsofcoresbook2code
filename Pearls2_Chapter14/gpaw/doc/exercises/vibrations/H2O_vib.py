"""Calculate the vibrational modes of a H2O molecule."""

from ase.vibrations import Vibrations
from gpaw import GPAW

h2o = GPAW('h2o.gpw', txt=None).get_atoms()

# Create vibration calculator
vib = Vibrations(h2o)
vib.run()
vib.summary(method='frederiksen')

# Make trajectory files to visualize normal modes:
for mode in range(9):
    vib.write_mode(mode)
