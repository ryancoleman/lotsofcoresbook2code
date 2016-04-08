# 12-ZGNR with adsorbed H at both edges.
# "ZGNR12.py" should be run first!

from ase.io import read
from ase import Atom
from ase.structure import graphene_nanoribbon
from ase.optimize import QuasiNewton
from gpaw import GPAW

GNR = read('ZGNR12.traj')
pos = GNR[22].position + [-1.05,0.0,0.0] 
GNR.append(Atom('H', pos))
pos = GNR[1].position + [1.05,0.0,0.0]
GNR.append(Atom('H', pos))
GNR.set_pbc((0,0,1))
kpts = (1,1,10)
calc = GPAW(kpts=kpts, spinpol=True)
GNR.set_calculator(calc)
dyn = QuasiNewton(GNR, trajectory='GNR_H.traj')
dyn.run(fmax=0.05)
