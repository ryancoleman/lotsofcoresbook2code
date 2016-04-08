from ase.optimize import BFGS
from ase.io import read, write
from ase.calculators.emt import EMT
from ase.ga.relax_attaches import VariansBreak
from ase.ga.atoms_attach import enable_raw_score_methods
import sys


fname = sys.argv[1]

print('Now relaxing {}'.format(fname))
a = read(fname)

a.set_calculator(EMT())
dyn = BFGS(a, trajectory=None, logfile=None)
vb = VariansBreak(a, dyn)
dyn.attach(vb.write)
dyn.run(fmax=0.05)

enable_raw_score_methods(a)
a.set_raw_score(-a.get_potential_energy())

write(fname[:-5] + '_done.traj', a)

print('Done relaxing {}'.format(fname))
