import ase
import ase.io as io
from ase.calculators.lj import LennardJones
from ase.constraints import FixBondLength
from ase.optimize import FIRE

for mic in [False, True]:
    a = ase.Atoms('CCC', positions=[[1,0,5],[0,1,5],[-1,0.5,5]], cell=[10,10,10], pbc=True)
    a.set_scaled_positions(a.get_scaled_positions()%1.0)
    a.set_calculator(LennardJones())
    a.set_constraint(FixBondLength(0, 2, mic=mic, atoms=a))

    dist = a.get_distance(0, 2, mic=mic)

    FIRE(a, logfile=None).run(fmax=0.01)

    assert abs(a.get_distance(0, 2, mic=mic) - dist) < 1e-6

