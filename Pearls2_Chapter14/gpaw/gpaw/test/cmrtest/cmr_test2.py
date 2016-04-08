# This test makes sure that the i/o interfaces work with CMR.
# CMR itself does not have to be installed for this test.
#
# The reason why CMR cannot use direct writes to DB/GPAW files is that 
# GPAW cannot always write a GPAW without performing a new calculation e.g. 
#    GPAW(filename).write(...)
# fails in some rare cases.

import os
from ase import Atom, Atoms
import gpaw.io
from gpaw import GPAW
from gpaw.test import equal
from gpaw.mpi import world, rank

import warnings
# cmr calls all available methods in ase.atoms detected by the module inspect.
# Therefore also deprecated methods are called - and we choose to silence those warnings.
warnings.filterwarnings('ignore', 'ase.atoms.*deprecated',)

a = 4.05
d = a / 2 ** 0.5
bulk = Atoms([Atom('Al', (0, 0, 0)),
              Atom('Al', (0.5, 0.5, 0.5))],
             pbc=True)
bulk.set_cell((d, d, a), scale_atoms=True)
h = 0.18
calc = GPAW(h=h,
            nbands=2 * 8,
            kpts=(2, 2, 2),
            convergence={'energy': 1e-5})
bulk.set_calculator(calc)
e0 = bulk.get_potential_energy()
calc.write("cmr_test2.gpw")

assert os.path.exists("cmr_test2.gpw")
reader = gpaw.io.open("cmr_test2.gpw", 'r')
w = {}
for key in reader.parameters:
    w[key] = reader.parameters[key]
for key in reader.shapes:
    w[key] = reader.get(key)
for key in reader.dims:
    w[key] = reader.dims[key]        
world.barrier()
if rank == 0:
    os.unlink("cmr_test2.gpw")
