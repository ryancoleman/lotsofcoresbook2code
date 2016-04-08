export MIC_ENV_PREFIX=MIC
export MIC_USE_2MB_BUFFERS=32K
export MIC_BUFFERSIZE=128M
export MIC_MKL_DYNAMIC=false
export MIC_OMP_NUM_THREADS=240
#export MIC_KMP_AFFINITY=verbose,granularity=fine,proclist=[1,2,3,4]
export MIC_KMP_AFFINITY=granularity=fine,balanced
