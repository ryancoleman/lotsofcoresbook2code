from __future__ import print_function
from ase import Atom
from gpaw import GPAW
from gpaw.cluster import Cluster
from gpaw.test import equal

h=0.40
txt = None
txt = 'mgga_sc.txt'

s = Cluster([Atom('H')])
s.minimal_box(4., h=h)
s.set_initial_magnetic_moments([1])
# see https://trac.fysik.dtu.dk/projects/gpaw/ticket/244
s.set_cell((8.400000,8.400000,8.400000))

c = GPAW(xc='TPSS', h=h, nbands=5, txt=txt, 
         eigensolver='rmm-diis',
         fixmom=True,
         maxiter=300)
c.calculate(s)

cpbe = GPAW(xc='PBE', h=h, nbands=5, txt=txt,
            eigensolver='rmm-diis',
            fixmom=True,
            maxiter=300)
cpbe.calculate(s)
cpbe.set(xc='TPSS')
cpbe.calculate()

print("Energy difference", (cpbe.get_potential_energy() - 
                            c.get_potential_energy()))
equal(cpbe.get_potential_energy(), c.get_potential_energy(), 0.002)
