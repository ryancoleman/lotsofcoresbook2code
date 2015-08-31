#!/bin/bash

# First, ENV variables:
# Edit or comment out these lines to suit your own installation:
export MIC_LD_LIBRARY_PATH=${MIC_LD_LIBRARY_PATH}:/ichec/packages/intel/intel-cluster-studio-2015/composer_xe_2015.0.090/compiler/lib/mic
# Add the MIC library to LD_LIBRARY_PATH to accomodate MIC_ENC_PREFIX: 
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/ichec/packages/intel/intel-cluster-studio-2015/composer_xe_2015.0.090/compiler/lib/mic

export I_MPI_PIN_DOMAIN=omp
export KMP_STACKSIZE=1g
export KMP_AFFINITY="scatter,granularity=fine"

export I_MPI_MIC=enable
export I_MPI_MIC_POSTFIX=.MIC

export MIC_ENV_PREFIX=MIC_
export KMP_STACKSIZE=200MB
export MIC_KMP_STACKSIZE=50MB
export MIC_KMP_MONITOR_STACKSIZE=12MB
export KMP_AFFINITY="scatter,granularity=fine"
export MIC_KMP_AFFINITY="scatter,granularity=fine"

# May need to cd to a directory containing executable & input data,
# ideally with "write" permission.
# cd <path-to-run-directory>

rm -f stresstest.dat

# ============

# First, a couple of runs on the Phi:
export OMP_NUM_THREADS=200
export MIC_OMP_NUM_THREADS=200
(time mpiexec.hydra -n 1 -hosts mic0 ./stest_hybrid_coarse.x < input_6GB_p1.txt ) >& t200_p1_6GB_phi.log 

export OMP_NUM_THREADS=100
export MIC_OMP_NUM_THREADS=100
(time mpiexec.hydra -n 2 -hosts mic0 ./stest_hybrid_coarse.x < input_6GB_p1x1x2.txt ) >& t100_p1x1x2_6GB_phi.log 

# ============

# Now a couple of runs on the host:
# 20 MPI processes on 20 physical cores, 2 threads/core (to use hyperthreading):
export OMP_NUM_THREADS=2
(time mpiexec.hydra -n 20 ./stest_hybrid_coarse.x < input_6GB_p2x2x5.txt ) >& t2_p2x2x5_6GB_host.log

# 40 MPI processes on all 40 logical cores (again using hyperthreads):
export OMP_NUM_THREADS=1
(time mpiexec.hydra -n 40 ./stest_hybrid_coarse.x < input_6GB_p2x4x5.txt ) >& t1_p2x4x5_6GB_host.log 

