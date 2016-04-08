
export MIC_ENV_PREFIX=MIC
export MIC_OMP_NUM_THREADS=240
export MIC_KMP_AFFINITY="explicit,proclist=[1-240],granularity=fine"
#export MIC_KMP_AFFINITY=scatter
export MIC_KMP_STACKSIZE=16m
export MIC_STACKSIZE=128M

