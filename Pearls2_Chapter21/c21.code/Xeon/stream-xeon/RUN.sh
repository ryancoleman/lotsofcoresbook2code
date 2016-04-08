#!/bin/bash

# 2 x Xeon E5-2697v2 12cores/CPU
#export OMP_NUM_THREADS=24
#export KMP_AFFINITY=proclist=[0-23:1],granularity=thread,explicit

# 2 x Xeon E5-2697v3 14cores/CPU
export OMP_NUM_THREADS=28
export KMP_AFFINITY=proclist=[0-27:1],granularity=thread,explicit

./stream_c.exe
