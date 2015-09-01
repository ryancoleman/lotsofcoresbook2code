#!/bin/bash
#
# ===== NOTE ========
# This script launches a timing_code on a Xeon
# 
# Either the OpenMP shared libs must be on the LD_LIBRARY_PATH automatically or they need
# to be set explicitly

# This is the hostname of my target Xeon Phi
# Your mileage may vary
TARGET_MIC=qcd12m

# This is the executable. In my case the host directory is mounted on the Xeon Phi
# Your mileage may vary
EXECUTABLE=`pwd`/build_avx_sse/tests/time_dslash
EXECUTABLE_ARGS=""

# Set the number of threads appropriate to your CPU
export OMP_NUM_THREADS=32
export KMP_AFFINITY=compact

# I use numactl to intereave memory allocations between my two sockets
# But feel free to use your favourite tool (likwid etc)
numactl --interleave=0,1 ${EXECUTABLE} ${EXECUTABLE_ARGS}
