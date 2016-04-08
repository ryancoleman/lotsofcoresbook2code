from ase.structure import nanotube
from gpaw import GPAW, Mixer, PoissonSolver, ConvergenceError
from gpaw.eigensolvers.rmm_diis import RMM_DIIS
from gpaw.mpi import size
from gpaw.mpi import rank
from gpaw.test import equal
from gpaw import use_mic

import time

tube = nanotube(6,6,10)

if rank == 0:
    print "Starting solver..."
tstart = time.time()
if use_mic:
    txt = 'out_nanotube_large_mic_p%d.txt' % size
else:
    txt = 'out_nanotube_large_p%d.txt' % size

conv = {'eigenstates' : 1e-4, 'density' : 1e-2, 'energy' : 1e-3}
calc = GPAW(# gpts=(96, 96, 128), 
            h=0.2,
            nbands=-60, width=0.1,
            poissonsolver=PoissonSolver(eps=1e-12),
            eigensolver=RMM_DIIS(keep_htpsit=False),
            # eigensolver=RMM_DIIS(),
            maxiter=6,
            mixer=Mixer(0.1, 5, 50),
            convergence=conv, txt=txt)
tube.set_calculator(calc)
try:
    e = tube.get_potential_energy()
except ConvergenceError:
    pass
tend = time.time()
if rank == 0:
    print "time for calculation: {0:.2f}sec".format(tend-tstart)
