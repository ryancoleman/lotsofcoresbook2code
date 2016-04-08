from __future__ import print_function
import numpy as np
from math import pi
from gpaw.coulomb import Coulomb
from gpaw.grid_descriptor import GridDescriptor
from gpaw.mpi import world, parallel
from gpaw.utilities.gauss import coordinates
from gpaw.test import equal
import time

def test_coulomb(N=2**6, a=20):
    Nc = (N, N, N)            # Number of grid point
    gd = GridDescriptor(Nc, (a, a, a), True)  # grid-descriptor object
    xyz, r2 = coordinates(gd) # matrix with the square of the radial coordinate
    r = np.sqrt(r2)           # matrix with the values of the radial coordinate
    nH = np.exp(-2 * r) / pi  # density of the hydrogen atom
    C = Coulomb(gd)           # coulomb calculator
    
    if parallel:
        C.load('real')
        t0 = time.time()
        print('Processor %s of %s: %s Ha in %s sec' % (
            gd.comm.rank + 1,
            gd.comm.size,
            -0.5 * C.coulomb(nH, method='real'),
            time.time() - t0))
        return
    else:
        C.load('recip_ewald')
        C.load('recip_gauss')
        C.load('real')
        test = {}
        t0 = time.time()
        test['dual density'] = (-0.5 * C.coulomb(nH, nH.copy()),
                                time.time() - t0)
        for method in ('real', 'recip_gauss', 'recip_ewald'):
            t0 = time.time()
            test[method] = (-0.5 * C.coulomb(nH, method=method),
                            time.time() - t0)
        return test

analytic = -5 / 16.0
res = test_coulomb(N=48, a=15)
if not parallel:
    print('Units: Bohr and Hartree')
    print('%12s %8s %8s' % ('Method', 'Energy', 'Time'))
    print('%12s %2.6f %6s' % ('analytic', analytic, '--'))
    for method, et in res.items():
        print('%12s %2.6f %1.7f' % ((method,) + et))

    equal(res['real'][0],         analytic, 6e-3)
    equal(res['recip_gauss'][0],  analytic, 6e-3)
    equal(res['recip_ewald'][0],  analytic, 2e-2)
    equal(res['dual density'][0], res['recip_gauss'][0], 1e-9)


# mpirun -np 2 python coulomb.py --gpaw-parallel --gpaw-debug
