import os
import numpy as np
from ase import Atom, Atoms
from gpaw import GPAW
from gpaw.test import equal
from gpaw.lrtddft import LrTDDFT

txt = '-'
txt = None
load = True
load = False
xc = 'LDA'

R = 0.7  # approx. experimental bond length
a = 4.0
c = 5.0
H2 = Atoms([Atom('H', (a / 2, a / 2, (c - R) / 2)),
            Atom('H', (a / 2, a / 2, (c + R) / 2))],
           cell=(a, a, c))

calc = GPAW(xc=xc, nbands=2, spinpol=False, eigensolver='rmm-diis', txt=txt)
H2.set_calculator(calc)
H2.get_potential_energy()
calc.write('H2saved_wfs.gpw', 'all')
calc.write('H2saved.gpw')
wfs_error = calc.wfs.eigensolver.error

#print "-> starting directly after a gs calculation"
lr = LrTDDFT(calc, txt='-')
lr.diagonalize()

#print "-> reading gs with wfs"
gs = GPAW('H2saved_wfs.gpw', txt=txt)

# check that the wfs error is read correctly, 
# but take rounding errors into account
assert( abs(calc.wfs.eigensolver.error/gs.wfs.eigensolver.error - 1) < 1e-8)
lr1 = LrTDDFT(gs, xc=xc, txt='-')
lr1.diagonalize()
# check the oscillator strrength
assert (abs(lr1[0].get_oscillator_strength()[0] /
            lr[0].get_oscillator_strength()[0] -1) < 1e-10)

#print "-> reading gs without wfs"
gs = GPAW('H2saved.gpw', txt=None)

lr2 = LrTDDFT(gs, txt='-')
lr2.diagonalize()
# check the oscillator strrength
d = abs(lr2[0].get_oscillator_strength()[0] /
        lr[0].get_oscillator_strength()[0] - 1)
assert (d < 2e-3), d
