#!/bin/bash 

mpiifort -openmp -g -traceback -O3 -mmic -o stest_hybrid_fine.x.MIC \
         stest_hybrid_fine.f  \
         shado_xchange_EW_fine.f  \
         shado_xchange_NS_fine.f  \
         shado_xchange_TB_fine.f

# These values are just placeholders, change as required:
export I_MPI_MIC=enable
export I_MPI_MIC_POSTFIX=.MIC
export I_MPI_PIN_DOMAIN=omp

export MIC_ENV_PREFIX=MIC_
export KMP_STACKSIZE=200MB
export MIC_KMP_STACKSIZE=50MB
export MIC_KMP_MONITOR_STACKSIZE=12MB
export KMP_AFFINITY="scatter,granularity=fine"
export MIC_KMP_AFFINITY="scatter,granularity=fine"

export OMP_NUM_THREADS=100
export MIC_OMP_NUM_THREADS=100
