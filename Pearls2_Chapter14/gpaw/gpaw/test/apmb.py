from __future__ import print_function
from ase import Atom, Atoms
from ase.parallel import parprint

from gpaw.test import equal
from gpaw import GPAW, mpi
from gpaw.lrtddft import LrTDDFT
import numpy

#numpy.seterr(all='raise')

txt='-'
txt='/dev/null'

load = False
#load = True
    
if not load:
    R=0.7 # approx. experimental bond length
    a = 3
    c = 4
    H2 = Atoms([Atom('H', (a/2,a/2,(c-R)/2)),
                Atom('H', (a/2,a/2,(c+R)/2))],
               cell=(a,a,c))
    calc = GPAW(xc='PBE', nbands=2, spinpol=False, txt=txt)
    H2.set_calculator(calc)
    H2.get_potential_energy()
##    calc.write('H2.gpw', 'all')
else:
    calc = GPAW('H2.gpw', txt=txt)
#calc.initialize_wave_functions()

#-----------------------------------------------------------
# DFT only

xc='LDA'

# no spin

lr = LrTDDFT(calc, xc=xc)
lr.diagonalize()

lr_ApmB = LrTDDFT(calc, xc=xc, force_ApmB=True)
lr_ApmB.diagonalize()
parprint('lr=', lr)
parprint('ApmB=', lr_ApmB)
equal(lr[0].get_energy(), lr_ApmB[0].get_energy(), 5.e-10)

# with spin
parprint('------ with spin')

if not load:
    c_spin = GPAW(xc='PBE', nbands=2, 
                  spinpol=True, parallel={'domain': mpi.world.size},
                  txt=txt)
    H2.set_calculator(c_spin)
    c_spin.calculate(H2)
##    c_spin.write('H2spin.gpw', 'all')
else:
    c_spin = GPAW('H2spin.gpw', txt=txt)
lr = LrTDDFT(c_spin, xc=xc)
lr.diagonalize()

lr_ApmB = LrTDDFT(c_spin, xc=xc, force_ApmB=True)
lr_ApmB.diagonalize()
parprint('lr=', lr)
parprint('ApmB=', lr_ApmB)
equal(lr[0].get_energy(), lr_ApmB[0].get_energy(), 5.e-10)
equal(lr[1].get_energy(), lr_ApmB[1].get_energy(), 5.e-10)

# with spin virtual
parprint('------ with virtual spin')

lr = LrTDDFT(calc, xc=xc, nspins=2)
lr.diagonalize()

# ApmB
lr_ApmB = LrTDDFT(calc, xc=xc, nspins=2)
lr_ApmB.diagonalize()
parprint('lr=', lr)
parprint('ApmB=', lr_ApmB)
equal(lr[0].get_energy(), lr_ApmB[0].get_energy(), 5.e-10)
equal(lr[1].get_energy(), lr_ApmB[1].get_energy(), 5.e-10)
    
#-------------------------------------------------------
# with HF exchange

xc='PBE0'

parprint('------ with spin xc=', xc)
lr_spin = LrTDDFT(c_spin, xc=xc)
lr_spin.diagonalize()
parprint('lr=', lr_spin)

parprint('------ with virtual spin xc=', xc)
lr = LrTDDFT(calc, xc=xc, nspins=2)
lr.diagonalize()
parprint('lr=', lr)
equal(lr[0].get_energy(), lr_spin[0].get_energy(), 3.8e-6)
equal(lr[1].get_energy(), lr_spin[1].get_energy(), 3.4e-6)
