"""Compare Al(fcc) and Al(bcc) at two different plane-wave cutoffs
and two differens k-point densities."""
from __future__ import print_function
from ase.lattice import bulk
from gpaw import GPAW, PW

afcc = 3.985
abcc = 3.190

for kdens in [2.0, 3.0]:
    for ecut in [300, 500]:
        fcc = bulk('Al', 'fcc', a=afcc)
        calc = GPAW(mode=PW(ecut),
                    kpts={'density': kdens},
                    txt='bulk-fcc-%.1f-%.1f.txt' % (ecut, kdens))
        fcc.set_calculator(calc)
        efcc = fcc.get_potential_energy()

        bcc = bulk('Al', 'bcc', a=abcc)
        calc = GPAW(mode=PW(ecut),
                    kpts={'density': 4.0},
                    txt='bulk-bcc-%.1f-%.1f.txt' % (ecut, kdens))
        bcc.set_calculator(calc)
        ebcc = bcc.get_potential_energy()

        print(kdens, ecut, efcc, ebcc, efcc - ebcc)
