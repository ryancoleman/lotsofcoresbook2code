from __future__ import print_function
from ase import Atom
from ase.neb import SingleCalculatorNEB
from ase.optimize import FIRE
from gpaw import GPAW
from gpaw.cluster import Cluster

txt=None

mol = Cluster([Atom('H'),
               Atom('H',[1,0,0]),
               Atom('H',[.5,.5,.5])],
              cell = [2,2,2],
              pbc=True)

def set_calculators(all=False):
    c=GPAW(h=.3, convergence={'eigenstates':0.1, 
                              'energy' : 0.1,
                              'density' : 0.01}, txt=txt)
#    c = EMT()
    n = len(images)
    if not all:
        n -= 2
    neb.set_calculators([c] * n)

images = [mol]
for i in range(4):
    images.append(images[0].copy())
images[-1].positions[2, 1] = 2 - images[0].positions[2, 1]
neb = SingleCalculatorNEB(images)
neb.interpolate()
for image in images:
    print(image[2].position)
set_calculators(True)

dyn = FIRE(neb, trajectory='mep.traj')
dyn.insert_observer(set_calculators)
print(dyn.run(fmax=8.))
