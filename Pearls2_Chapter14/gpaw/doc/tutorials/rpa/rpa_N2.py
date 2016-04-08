from __future__ import print_function
from ase.parallel import paropen
from ase.units import Hartree
from gpaw.xc.rpa import RPACorrelation

rpa1 = RPACorrelation('N.gpw', nblocks=8, wstc=True, txt='rpa_N.txt')
rpa2 = RPACorrelation('N2.gpw', nblocks=8, wstc=True, txt='rpa_N2.txt')

E1_i = rpa1.calculate(ecut=400)
E2_i = rpa2.calculate(ecut=400)

f = paropen('rpa_N2.dat', 'w')
for ecut, E1, E2 in zip(rpa1.ecut_i, E1_i, E2_i):
    print(ecut * Hartree, E2 - 2 * E1, file=f)
f.close()
