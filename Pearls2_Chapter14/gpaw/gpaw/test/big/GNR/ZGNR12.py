from ase.optimize import QuasiNewton
from ase.structure import graphene_nanoribbon
from gpaw import GPAW

GNR = graphene_nanoribbon(12,1,type='zigzag', vacc=6)
GNR.set_pbc((0,0,1))
kpts = (1,1,10)
calc = GPAW(kpts=kpts,spinpol=True)
GNR.set_calculator(calc)
dyn = QuasiNewton(GNR, trajectory='ZGNR12.traj')
dyn.run(fmax=0.05)
