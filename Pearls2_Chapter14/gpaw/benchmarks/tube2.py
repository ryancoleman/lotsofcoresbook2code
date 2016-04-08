from ase.structure import nanotube
from gpaw import GPAW, Mixer, PoissonSolver
from gpaw.mpi import size
from gpaw.test import equal
from gpaw import use_mic
from gpaw.mpi import rank

import time

tube = nanotube(6,6,9)

if rank == 0:
    print "Starting solver..."
tstart = time.time()
if use_mic:
    txt = 'out_nanotube2_mic_p%d.txt' % size
else:
    txt = 'out_nanotube2_p%d.txt' % size
conv = {'eigenstates' : 1e-4, 'density' : 1e-2, 'energy' : 1e-3}
calc = GPAW(gpts=(96, 96, 128), nbands=-20, width=0.1,
            poissonsolver=PoissonSolver(eps=1e-12),
            mixer=Mixer(0.1, 5, 50),
            convergence=conv, txt=txt)
tube.set_calculator(calc)
e = tube.get_potential_energy()
tend = time.time()
if rank == 0:
    print "time for calculation: {0:.2f}sec".format(tend-tstart)
equal(e, -2159.016150, 1e-2)
