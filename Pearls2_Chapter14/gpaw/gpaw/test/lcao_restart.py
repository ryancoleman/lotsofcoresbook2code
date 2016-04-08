from __future__ import print_function
# Test that LCAO wavefunctions are available and equal after restarts
# in normal as well as 'all' mode

import numpy as np
from ase.structure import molecule
from gpaw import GPAW

# setting number of decimals globally makes numpy.test() tests
# which use docstrings fail
#np.set_printoptions(precision=3, suppress=1)

system = molecule('H2')
system.center(vacuum=2.5)

calc = GPAW(mode='lcao', basis='sz(dzp)', h=0.3, nbands=1, txt=None)
system.set_calculator(calc)
system.get_potential_energy()
wf = calc.get_pseudo_wave_function(0)

for mode in ['normal', 'all']:
    fname = 'lcao-restart.%s.gpw' % mode
    calc.write(fname, mode=dict(normal='', all='all')[mode])
    
    calc2 = GPAW(fname, txt=None)
    wf2 = calc2.get_pseudo_wave_function(0)
    err = np.abs(wf2 - wf).max()
    print('%s: err=%s' % (mode, repr(err)))
    assert err == 0.0
