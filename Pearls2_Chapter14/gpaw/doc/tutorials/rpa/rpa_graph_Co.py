from __future__ import print_function
from ase.parallel import paropen
from gpaw.xc.rpa import RPACorrelation

ds = [1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 5.0, 6.0, 10.0]

for d in ds:
    rpa = RPACorrelation('gs_%s.gpw' % d, txt='rpa_%s_%s.txt' % (ecut, d))
    E_rpa = rpa.calculate(ecut=[200],
                          frequency_scale=2.5,
                          skip_gamma=True,
                          filename='restart_%s_%s.txt' % (ecut, d))

    f = paropen('rpa_%s.dat' % ecut, 'a')
    print(d, E_rpa, file=f)
    f.close()
