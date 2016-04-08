from __future__ import print_function
import sys

from ase import Atoms
from ase.utils.timing import Timer

from gpaw import GPAW
from gpaw.test import equal
from gpaw.xc.hybrid import HybridXC

timer = Timer()

loa = Atoms('Be2',
            [(0, 0, 0), (2.45, 0, 0)],
            magmoms=[0.5, 0.5],
            cell=[5.9, 4.8, 5.0])
loa.center()

fgl = [False, True]
#fgl = [True, False]

txt='-'
txt='/dev/null'

E = {}
niter = {}
for fg in fgl:
    if fg:
        tstr = 'Exx on fine grid'
    else:
        tstr = 'Exx on coarse grid'
    timer.start(tstr)
    calc = GPAW(h=0.3,
                eigensolver='rmm-diis',
                xc='PBE',
                nbands=4,
                convergence={'eigenstates': 1e-4},
                charge=-1)
    loa.set_calculator(calc)
    E[fg] = loa.get_potential_energy()
    calc.set(xc=HybridXC('PBE0', finegrid=fg))
    E[fg] = loa.get_potential_energy()
    niter[fg] = calc.get_number_of_iterations()
    timer.stop(tstr)

timer.write(sys.stdout)

print('Total energy on the fine grid   =', E[True])
print('Total energy on the coarse grid =', E[False])
equal(E[True], E[False], 0.01)

energy_tolerance = 0.0003
equal(E[False], 6.97818, energy_tolerance)
equal(E[True], 6.97153, energy_tolerance)
