# creates: pt_h2.png

from ase import Atoms
from ase.io import write

a = 2.41  # Pt binding lenght
b = 0.90  # H2 binding lenght
c = 1.70  # Pt-H binding lenght
L = 7.00  # width of unit cell

# Setup the Atoms for the scattering region.
atoms = Atoms('Pt5H2Pt5',
              cell=[3 * a + b + 2 * c, L, L],
              pbc=True)

atoms.positions[:5, 0] = [(i - 2.5) * a for i in range(5)]
atoms.positions[-5:, 0] = [(i - 2.5) * a + b + 2 * c for i in range(4, 9)]
atoms.positions[5:7, 0] = [1.5 * a + c, 1.5 * a + c + b]
atoms.positions[:, 1:] = L / 2.

write('pt_h2.pov', atoms, show_unit_cell=2,
      transparent=False, display=False, run_povray=True)
