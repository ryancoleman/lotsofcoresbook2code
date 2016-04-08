from ase.structure import nanotube
from gpaw import GPAW, RMM_DIIS
from gpaw.mpi import size
from gpaw.test import equal
from gpaw import use_mic
from gpaw.mpi import rank

import time

tube = nanotube(8,8,8)

if rank == 0:
    print "Starting solver..."
tstart = time.time()
if use_mic:
    txt = 'out_nanotube1_mic_p%d.txt' % size
else:     
    txt = 'out_nanotube1_p%d.txt' % size
conv = {'eigenstates' : 1e-4, 'density' : 1e-2, 'energy' : 1e-2}
calc = GPAW(h=0.2, nbands=-20, width=0.1, kpts=(1,1,2), 
            eigensolver=RMM_DIIS(keep_htpsit=False),
            #convergence=conv, txt=txt)
            convergence=conv)
tube.set_calculator(calc)
e = tube.get_potential_energy()
tend = time.time()
if rank == 0:
    print "time for calculation: {0:.2f}sec".format(tend-tstart)
equal(e, -628.838671, 1e-2)
