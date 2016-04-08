from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from gpaw.wavefunctions.pw import PW


# Plane wave cutoff
pwcutoff = 400.0

# NxNxN k-point sampling
k = 4

# Si lattice constant
alat = 5.421

# bulk calculation
bulk_crystal = bulk('Si', 'diamond', a=alat)
bulk_calc = GPAW(mode=PW(pwcutoff),
                 kpts={'size': (k, k, k), 'gamma': True},  # gamma-centred grid
                 xc='PBE',
                 occupations=FermiDirac(0.01),
                 txt='si.rpa.pbe_output.txt',
                 parallel={'band': 1})

bulk_crystal.set_calculator(bulk_calc)
e0_bulk_pbe = bulk_crystal.get_potential_energy()

# Now we have the density, but only the occupied states;
# at this point we diagonalise the full Hamiltonian to get
# the empty states as well (expensive)
bulk_calc.diagonalize_full_hamiltonian(nbands=200)

# the 'all' ensures we write wavefunctions too
bulk_calc.write('bulk.all.gpw', mode='all')
