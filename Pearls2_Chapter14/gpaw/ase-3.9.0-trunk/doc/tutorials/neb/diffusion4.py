from ase.io import read
from ase.calculators.emt import EMT
from ase.neb import NEB
from ase.optimize import BFGS

# read the last structures (of 5 images used in NEB)
images = read('neb.traj@-5:')

for i in range(1, len(images) - 1):
    images[i].set_calculator(EMT())

neb = NEB(images)
qn = BFGS(neb, trajectory='neb_restart.traj')
qn.run(fmax=0.005)
