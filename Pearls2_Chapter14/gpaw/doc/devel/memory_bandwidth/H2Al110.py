#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser

code_choices = ['gpaw', 'dacapo']

parser = OptionParser(usage='%prog [options] package.\nExample of call:\n'+
                      'python %prog\n',
                      version='%prog 0.1')
parser.add_option('--code', dest="code", type="choice",
                  default=code_choices[0],
                  choices=code_choices,
                  help='code: which code to use.')
parser.add_option("--runs", dest="runs",
                  default=7,
                  help='use that many runs to calculate the average.')
parser.add_option('-v', '--verbose', action='store_true',
                  default=False,
                  help='verbose mode.')


opt, args = parser.parse_args()

from os import remove
from os.path import exists

try:
    import numpy as np
except ImportError:
    raise SystemExit('numpy is not installed!')

try:
    import gpaw
except ImportError:
    raise SystemExit('gpaw is not installed!')

from gpaw.utilities.tools import gridspacing2cutoff

try:
    import ase
except ImportError:
    raise SystemExit('ase is not installed!')

from ase import Atoms, Atom

import time

a = 4.00
d = a / 2**0.5
z = 1.1
b = 1.5

def memory_bandwidth(code='gpaw', runs=7):

    slab = Atoms([Atom('Al', (0, 0, 0)),
                        Atom('Al', (a, 0, 0)),
                        Atom('Al', (a/2, d/2, -d/2)),
                        Atom('Al', (3*a/2, d/2, -d/2)),
                        Atom('Al', (0, 0, -d)),
                        Atom('Al', (a, 0, -d)),
                        Atom('Al', (a/2, d/2, -3*d/2)),
                        Atom('Al', (3*a/2, d/2, -3*d/2)),
                        Atom('Al', (0, 0, -2*d)),
                        Atom('Al', (a, 0, -2*d)),
                        Atom('H', (a/2-b/2, 0, z)),
                        Atom('H', (a/2+b/2, 0, z))],
                       cell=(2*a, d, 5*d), pbc=(1, 1, 1))

    h = 0.15
    nbands = 28
    kpts = (2, 6, 1)

    parameters = {}

    if code == 'gpaw':
        from gpaw import Calculator
        from gpaw.mpi import rank
        parameters['convergence'] = {'eigenstates': 1e-5}
        parameters['h'] = h
    elif code == 'dacapo':
        from ase.calculators.dacapo import Dacapo as Calculator
        parameters['planewavecutoff'] = gridspacing2cutoff(h)
        parameters['densitycutoff'] = parameters['planewavecutoff']*1.5
        rank = 0

    t = 0.0
    t_runs = []
    for n in range(runs):
        t0 = time.time()
        for i in range(1):
            calc = Calculator(
                nbands=nbands,
                kpts=kpts,
                **parameters)
            slab.set_calculator(calc)
            e = slab.get_potential_energy()
            del calc
            if exists('out.nc'): remove('out.nc')
        t1 = time.time()
        t = t + t1 - t0
        t_runs.append(t1 - t0)
        print('Run: ', n, ' energy ', e, ' rank: ', str(rank), ' time: ', time.time() - t0)
    if rank == 0:
        print('Rank '+str(rank)+': time [sec]: avg '+str(round(np.average(t_runs),1))+', stddev '+str(round(np.std(t_runs),1))+', min '+str(round(min(t_runs),1))+', max '+str(round(max(t_runs),1)))


if __name__ == '__main__':


    code = opt.code
    assert code in code_choices, code+' not in '+str(code_choices)

    if code == 'dacapo':
        try:
            import ASE
        except ImportError:
            raise SystemExit('ASE (2) is not installed!')

    runs = int(opt.runs)
    assert runs >= 1, runs+' must be >= 1'

    memory_bandwidth(code=code, runs=runs)
