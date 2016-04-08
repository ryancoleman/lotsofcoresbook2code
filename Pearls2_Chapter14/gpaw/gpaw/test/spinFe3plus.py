from ase import Atoms, Atom
import numpy as np

from gpaw import GPAW, FermiDirac
from gpaw.test import equal

h=.4
q=3
spin=True

s = Atoms([Atom('Fe')])
s.center(vacuum=2.5)
convergence={'eigenstates':0.01, 'density':0.1, 'energy':0.1}

# use Hunds rules

c = GPAW(charge=q, h=h, nbands=5,
         hund=True,
         eigensolver='rmm-diis',
         occupations=FermiDirac(width=0.1),
         convergence=convergence
         )
c.calculate(s)
equal(c.get_magnetic_moment(0), 5, 0.1)

# set magnetic moment

mm = [5]
s.set_initial_magnetic_moments(mm)
c = GPAW(charge=q, h=h, nbands=5,
         occupations=FermiDirac(width=0.1, fixmagmom=True),
         convergence=convergence
         )
c.calculate(s)
equal(c.get_magnetic_moment(0), 5, 0.1)
