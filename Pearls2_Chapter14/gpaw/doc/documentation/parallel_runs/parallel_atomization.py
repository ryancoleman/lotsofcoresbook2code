"""This script calculates the atomization energy of nitrogen using two
processes, each process working on a separate system."""
from __future__ import print_function
from gpaw import GPAW, mpi
import numpy as np
from ase import Atoms, Atom

cell = (8., 8., 8.)
p = 4.
separation = 1.103

rank = mpi.world.rank

# Master process calculates energy of N, while the other one takes N2
if rank == 0:
    system = Atoms('N', [(p, p, p)], magmoms=[3], cell=cell)
elif rank == 1:
    system = Atoms('N2', [(p, p, p + separation / 2.),
                          (p, p, p - separation / 2.)], 
                   cell=cell)
else:
    raise Exception('This example uses only two processes')

# Open different files depending on rank
output = '%d.txt' % rank
calc = GPAW(communicator=[rank], txt=output, xc='PBE')
system.set_calculator(calc)
energy = system.get_potential_energy()

# Now send the energy from the second process to the first process,
if rank == 1:
    # Communicators work with arrays from Numeric only:
    mpi.world.send(np.array([energy]), 0)
else:
    # The first process receives the number and prints the atomization energy
    container = np.array([0.])
    mpi.world.receive(container, 1)

    # Ea = E[molecule] - 2 * E[atom]
    atomization_energy = container[0] - 2 * energy
    print('Atomization energy: %.4f eV' % atomization_energy)
