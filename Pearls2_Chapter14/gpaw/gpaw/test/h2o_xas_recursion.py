from __future__ import print_function
import os
from math import pi, cos, sin
from ase import Atom, Atoms
from ase.parallel import rank, barrier
from gpaw import GPAW
from gpaw.test import equal, gen
import numpy as np

# Generate setup for oxygen with half a core-hole:
gen('O', name='hch1s', corehole=(1, 0, 0.5))

if 1:
    a = 5.0
    d = 0.9575
    t = pi / 180 * 104.51
    H2O = Atoms([Atom('O', (0, 0, 0)),
                 Atom('H', (d, 0, 0)),
                 Atom('H', (d * cos(t), d * sin(t), 0))],
                cell=(a, a, a), pbc=False)
    H2O.center()
    calc = GPAW(nbands=10, h=0.2, setups={'O': 'hch1s'})
    H2O.set_calculator(calc)
    e = H2O.get_potential_energy()
    niter = calc.get_number_of_iterations()
    calc.write('h2o.gpw')
else:
    calc = GPAW('h2o.gpw')
    calc.initialize_positions()

from gpaw.xas import RecursionMethod

if 1:
    r = RecursionMethod(calc)
    r.run(400)
    r.write('h2o.pckl')
else:
    r = RecursionMethod(filename='h2o.pckl')

if 0:
    from pylab import *

    x = -30 + 40 * np.arange(300) / 300.0

    for n in range(50, 401, 50):
        y = r.get_spectra(x, imax=n)
        plot(x, y[0], label=str(n))
    legend()
    show()

print(e, niter)
energy_tolerance = 0.0002
niter_tolerance = 0
equal(e, -17.9621, energy_tolerance)
