from ase import Atoms
from gpaw import GPAW

# Sodium dimer, Na2
d = 1.5
atoms = Atoms(symbols='Na2',
              positions=[(0, 0, d),
                         (0, 0, -d)],
              pbc=False)

# For real calculations use larger vacuum (e.g. 6)
atoms.center(vacuum=4.0)

calc = GPAW(nbands=1,
            h=0.35,
            setups={'Na': '1'},
            txt='Na2_gs.txt')

atoms.set_calculator(calc)
e = atoms.get_potential_energy()

# Calculate also unoccupied states with the fixed density
# unoccupied states converge often better with cg
calc.set(nbands=20, convergence={'bands': 20},
         eigensolver='cg',
         fixdensity=True)
e = atoms.get_potential_energy()
# write the wave functions to a file
calc.write('na2_gs_unocc.gpw', 'all')
