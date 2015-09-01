#!/bin/bash

unset MIC_OMP_PROC_BIND
unset MIC_OMP_PLACES

export OFFLOAD_DEVICES=0 
export OFFLOAD_INIT=on_offload 
export MIC_ENV_PREFIX=MIC 
export MIC_KMP_AFFINITY=compact,granularity=fine
export MIC_OMP_NESTED=FALSE

mpirun -np 1 part-2/modal_man parameter_DDX9_smica_2000_4_7_601.ini
