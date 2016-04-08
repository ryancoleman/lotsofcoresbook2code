from __future__ import print_function
import sys
import numpy as np

from ase.structure import molecule
from ase.parallel import parprint

from gpaw import GPAW, PW
from gpaw.analyse.multipole import Multipole
from gpaw.cluster import Cluster
from gpaw.test import equal

h = 0.3

s = Cluster(molecule('H2O'))
s.minimal_box(3., h)

gpwname = 'H2O_h' + str(h) + '.gpw'
try:
    # XXX check why this fails in parallel
    calc = GPAW(gpwname + 'failsinparallel', txt=None)
    atoms = calc.get_atoms()
    calc.density.set_positions(atoms.get_scaled_positions() % 1.0)
    calc.density.interpolate_pseudo_density()
    calc.density.calculate_pseudo_charge()
except IOError:
    calc = GPAW(h=h, charge=0, txt=None)
    calc.calculate(s)
    calc.write(gpwname)

dipole_c = calc.get_dipole_moment()
parprint('Dipole', dipole_c) 

center = np.array([1,1,1]) * 50.
mp = Multipole(center, calc, lmax=2)
q_L = mp.expand(-calc.density.rhot_g)
parprint('Multipole', q_L)

# The dipole moment is independent of the center
equal(dipole_c[2], q_L[2], 1e-10)

mp.to_file(calc, mode='w')
