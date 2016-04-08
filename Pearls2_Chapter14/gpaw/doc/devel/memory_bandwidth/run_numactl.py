import os

from gpaw.mpi import rank
machine = os.environ.get('MACHINE', 'TEST')
ncores = os.environ.get('NCORES', 8)
if rank == 0:
    os.chdir(machine)
    os.system('STARTCORES=%d &&. ../run_numactl.sh' % ncores)
    #os.system('. ../run_numactl.sh') # full memory benchmark
