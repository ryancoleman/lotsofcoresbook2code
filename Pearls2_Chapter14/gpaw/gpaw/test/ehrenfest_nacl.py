from __future__ import print_function
from ase import Atoms
from gpaw import GPAW
from gpaw.tddft import TDDFT
from gpaw.tddft.ehrenfest import EhrenfestVelocityVerlet
from gpaw.test import equal

d = 4.0
atoms = Atoms('NaCl', [(0,0,0),(0,0,d)])
atoms.center(vacuum=4.5)

gs_calc = GPAW(nbands=4, eigensolver='cg', gpts=(32, 32, 44), xc='LDA',
               setups={'Na': '1'})
atoms.set_calculator(gs_calc)
atoms.get_potential_energy()

gs_calc.write('nacl_gs.gpw', 'all')

td_calc = TDDFT('nacl_gs.gpw', propagator='EFSICN')
evv = EhrenfestVelocityVerlet(td_calc, 0.001)

i=0
evv.get_energy()
r = evv.x[1][2] - evv.x[0][2]
# print 'E = ', [i, r, evv.Etot, evv.Ekin, evv.Epot]
    
for i in range(5):
    evv.propagate(1.0)
    evv.get_energy()
    r = evv.x[1][2] - evv.x[0][2]
    print('E = ', [i+1, r, evv.Etot, evv.Ekin, evv.Epot])

equal(r, 7.558883144, 1e-7)
equal(evv.Etot, -0.1036763317, 1e-7)
