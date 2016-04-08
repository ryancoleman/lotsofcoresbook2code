import os
from ase.vibrations import Vibrations
from gpaw import GPAW
from gpaw.mpi import world

h2o = GPAW('h2o.gpw', txt=None).get_atoms()

# Test restart:
vib = Vibrations(h2o)
vib.run()
vib.summary(method='frederiksen')

# Remove a displacement file and redo it:
if world.rank == 0:
    os.remove('vib.1z-.pckl')
world.barrier()
vib = Vibrations(h2o)
vib.run()
vib.summary(method='frederiksen')
