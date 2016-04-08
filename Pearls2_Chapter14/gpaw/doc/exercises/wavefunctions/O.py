from ase import Atoms
from ase.io import write
from gpaw import GPAW

# Oxygen atom:
atom = Atoms('O', cell=[6, 6, 6], pbc=False)
atom.center()

calc = GPAW(h=0.2,
            hund=True,  # assigns the atom its correct magnetic moment
            txt='O.txt')

atom.set_calculator(calc)
atom.get_potential_energy()

# Write wave functions to gpw file:
calc.write('O.gpw', mode='all')

# Generate cube-files of the orbitals:
for spin in [0, 1]:
    for n in range(calc.get_number_of_bands()):
        wf = calc.get_pseudo_wave_function(band=n, spin=spin)
        write('O.%d.%d.cube' % (spin, n), atom, data=wf)
