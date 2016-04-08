import os
import shutil

from gpaw.mpi import rank
machine = os.environ.get('MACHINE', 'TEST')
ncores = os.environ.get('NCORES', 8)
if rank == 0:
    os.chdir(machine)
    os.system('python ../memory_bandwidth.py --runs=5 --startcores='+str(ncores))
    #os.system('python ../memory_bandwidth.py --runs=5') # full memory benchmark
    shutil.copy('memory_bandwidth_'+machine+'_py.png', '..')
