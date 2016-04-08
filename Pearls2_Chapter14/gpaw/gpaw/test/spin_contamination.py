import os
import sys

from gpaw.test import equal
from ase import Atom

from gpaw import GPAW, MixerDif
from gpaw.cluster import Cluster

h = .25
q = 0
box = 3.
spin=True

# B should not have spin contamination
s = Cluster([Atom('B')])
s.minimal_box(box, h=h)
s.set_initial_magnetic_moments([-1])

c = GPAW(xc='LDA', nbands=-3, 
         charge=q, spinpol=spin, h=h,
         mixer=MixerDif(beta=0.05, nmaxold=5, weight=50.0),
         convergence={'eigenstates': 0.078, 'density': 1e-2, 'energy': 0.1},
         )
c.calculate(s)
equal(c.density.get_spin_contamination(s, 1), 0., 0.01) 

# setup H2 at large distance with different spins for the atoms
s = Cluster([Atom('H'), Atom('H',[0,0,3.0])])
s.minimal_box(box, h=h)
s.set_initial_magnetic_moments([-1,1])

c = GPAW(xc='LDA', nbands=-3, 
         charge=q, spinpol=spin, h=h,
         convergence={'eigenstates': 0.078, 'density': 1e-2, 'energy': 0.1},
         )
c.calculate(s)

scont_s = [c.density.get_spin_contamination(s), 
           c.density.get_spin_contamination(s, 1)]
equal(scont_s[0], scont_s[1], 2.e-4) # symmetry
equal(scont_s[0], 0.9655, 1.5e-3)
