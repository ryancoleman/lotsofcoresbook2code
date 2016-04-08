from __future__ import print_function
from ase.parallel import paropen
from ase.units import Hartree
from gpaw.xc.rpa import RPACorrelation

f = paropen('con_freq.dat', 'w')
for N in [4, 6, 8, 12, 16, 24, 32]:
    rpa = RPACorrelation('N2.gpw', txt='rpa_N2_frequencies.txt', nfrequencies=N)
    E = rpa.calculate(ecut=[50])
    print(N, E[0], file=f)
    if N == 16:
        f16 = paropen('frequency_gauss16.dat', 'w')
        for w, e in zip(rpa.omega_w, rpa.E_w):
            print(w * Hartree, e, file=f16)
        f16.close()
f.close()
