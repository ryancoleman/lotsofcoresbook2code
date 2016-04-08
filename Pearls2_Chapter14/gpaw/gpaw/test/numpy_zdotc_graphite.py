# gpaw-python Segmentation faults
# when gpaw-python and numpy are linked to different blas

import numpy as np
import sys

from math import sqrt
from ase import Atom, Atoms
from gpaw import GPAW
from gpaw import ConvergenceError

kpts = (2,1,1)
a=1.42
c=3.355

# AB stack
atoms = Atoms('C4',[
              (1/3.0,1/3.0,0),
              (2/3.0,2/3.0,0),
              (0.   ,0.   ,0.5),
              (1/3.0,1/3.0,0.5)
              ],
              pbc=(1,1,1))
atoms.set_cell([(sqrt(3)*a/2.0,3/2.0*a,0),
                (-sqrt(3)*a/2.0,3/2.0*a,0),
                (0.,0.,2*c)],
               scale_atoms=True)

calc = GPAW(gpts=(8, 8, 20), nbands=9, kpts=kpts, maxiter=1)

atoms.set_calculator(calc)

try:
    pot = atoms.get_potential_energy()
except ConvergenceError:
    pass
