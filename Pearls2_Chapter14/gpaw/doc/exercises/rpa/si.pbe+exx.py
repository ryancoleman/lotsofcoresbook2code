from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from gpaw.wavefunctions.pw import PW
from gpaw.xc.exx import EXX

from ase.parallel import paropen

resultfile = paropen('si.pbe+exx.results.txt', 'a')

# Plane wave cutoff
pwcutoff = 400.0

# NxNxN k-point sampling, gamma-centred grid
k = 6

# Si lattice constant
alat = 5.421

# Do the bulk calculation

bulk_crystal = bulk('Si', 'diamond', a=alat)
bulk_calc = GPAW(mode=PW(pwcutoff),
                 kpts={'size': (k, k, k), 'gamma': True},
                 xc='PBE',
                 occupations=FermiDirac(0.01),
                 txt='si.pbe+exx.pbe_output.txt',
                 parallel={'band': 1}
                 )

bulk_crystal.set_calculator(bulk_calc)
e0_bulk_pbe = bulk_crystal.get_potential_energy()

#  Write to file
bulk_calc.write('bulk.gpw', mode='all')

# Now the exact exchange
exx_bulk = EXX('bulk.gpw', txt='si.pbe+exx.exx_output.txt')
exx_bulk.calculate()
e0_bulk_exx = exx_bulk.get_total_energy()

s = str(alat)
s += ' '
s += str(k)
s += ' '
s += str(pwcutoff)
s += ' '
s += str(e0_bulk_pbe)
s += ' '
s += str(e0_bulk_exx)
s += '\n'
resultfile.write(s)
