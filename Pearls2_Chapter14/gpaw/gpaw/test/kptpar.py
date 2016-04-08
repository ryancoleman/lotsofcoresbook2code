# see changeset 4891
import numpy as np
from ase import Atoms
from gpaw import GPAW
from gpaw.mpi import world
from gpaw.test import equal

a = 2.5
H = Atoms('H', cell=[a, a, a], pbc=True)

energy_tolerance = 0.00006
niter_tolerance = 0

if world.size >= 3:
    calc = GPAW(kpts=[6, 6, 1],
                spinpol=True,
                parallel={'domain': world.size},
                txt='H-a.txt')
    H.set_calculator(calc)
    e1 = H.get_potential_energy()
    niter1 = calc.get_number_of_iterations()
    assert H.get_calculator().wfs.kd.comm.size == 1

    equal(e1, -2.23708481, energy_tolerance)

    if world.rank < 3:
        comm = world.new_communicator(np.array([0, 1, 2]))
        H.set_calculator(GPAW(kpts=[6, 6, 1],
                              spinpol=True,
                              communicator=comm,
                              txt='H-b.txt'))
        e2 = H.get_potential_energy()
        assert H.get_calculator().wfs.kd.comm.size == 3
        equal(e1, e2, 5e-9)
    else:
        comm = world.new_communicator(np.array(range(3, world.size)))
        H.set_calculator(GPAW(kpts=[6, 6, 1],
                              spinpol=True,
                              communicator=comm,
                              txt='H-b2.txt'))
        e2 = H.get_potential_energy()
