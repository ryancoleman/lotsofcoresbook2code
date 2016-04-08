from ase.lattice.spacegroup import crystal
from gpaw import GPAW, PW

a = 4.0351  # Experimental lattice constant in Angstrom
Ecut = 250  # Energy cut off for PW calculation
k = 4  # Number of kpoints per each direction

# This gives the typical NaCl structure:
LiF = crystal(['Li', 'F'],
              [(0, 0, 0), (0.5, 0.5, 0.5)],
              spacegroup=225,
              cellpar=[a, a, a, 90, 90, 90])

calc = GPAW(mode=PW(Ecut),
            xc='LDA',
            kpts=(k, k, k),
            txt='LiF_out_gs.txt')

LiF.set_calculator(calc)
LiF.get_potential_energy()

# With the full diagonalization we calculate all the single particle states
# needed for the response calculation:
calc.diagonalize_full_hamiltonian(nbands=100)
calc.write('LiF_fulldiag.gpw', 'all')
