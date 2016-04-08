from __future__ import print_function
import numpy as np
from ase import Atoms
from gpaw import GPAW
from gpaw.test import equal

bulk = Atoms('Li', pbc=True)
k = 4
g = 8
calc = GPAW(gpts=(g, g, g), kpts=(k, k, k),
                  mode='lcao', basis='dzp')
bulk.set_calculator(calc)
e = []
niter = []
A = [2.6, 2.65, 2.7, 2.75, 2.8]
for a in A:
    bulk.set_cell((a, a, a))
    e.append(bulk.get_potential_energy())
    niter.append(calc.get_number_of_iterations())

a = np.roots(np.polyder(np.polyfit(A, e, 2), 1))[0]
print('a =', a)
equal(a, 2.63781, 0.0001)

e_ref = [-1.8677343236247692, -1.8690343169380492, -1.8654175796625045,
         -1.8566274574918875, -1.8432374955346396]
niter_ref = [6, 6, 6, 6, 6]

print(e)
energy_tolerance = 0.00003
niter_tolerance = 0

for i in range(len(A)):
    equal(e[i], e_ref[i], energy_tolerance)
    equal(niter[i], niter_ref[i], niter_tolerance)

wf1 = calc.get_pseudo_wave_function(kpt=3, band=0)
calc.write('Li', mode='all')
calc2 = GPAW('Li')
calc2.initialize_positions()
wf2 = calc2.get_pseudo_wave_function(kpt=3, band=0)
equal(abs(wf1 - wf2).max(), 0, 1e-9)

