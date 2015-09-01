#!/bin/bash

export OFFLOAD_DEVICES=0 
export OFFLOAD_INIT=on_offload 
export MIC_ENV_PREFIX=MIC 
export MIC_KMP_AFFINITY= 
export MIC_OMP_PROC_BIND="spread,close" 
export MIC_OMP_PLACES=threads 
export MIC_OMP_NESTED=FALSE

# experiment with values {default, 0, infinite}
export MIC_KMP_BLOCKTIME=infinite

mpirun -np 1 part-2/modal_teams parameter_DDX9_smica_2000_4_7_601.ini
