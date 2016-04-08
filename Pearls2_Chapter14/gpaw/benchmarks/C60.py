from ase.data.extra_molecules import data
from ase.structure import molecule
from gpaw import GPAW, ConvergenceError
from gpaw.eigensolvers import RMM_DIIS
from gpaw.mpi import size
from _gpaw import offload_report
from gpaw import use_mic
from gpaw.mpi import rank

import time

offload_report(0)

if use_mic:
    txt = 'out_C60_mic_p%d.txt' % size
else:
    txt = 'out_C60_p%d.txt' % size

if rank == 0:
    print "Starting solver..."
tstart = time.time()
atoms = molecule('C60', data=data)
atoms.center(3.5)
calc = GPAW(h=0.18, nbands=400, eigensolver=RMM_DIIS(keep_htpsit=False),
            txt=txt,
            maxiter=10)
atoms.set_calculator(calc)
try:
    atoms.get_potential_energy()
except ConvergenceError:
    pass
tend = time.time()
if rank == 0:
    print "time for calculation: {0:.2f}sec".format(tend-tstart)
