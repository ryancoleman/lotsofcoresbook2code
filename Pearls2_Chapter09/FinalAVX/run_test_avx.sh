#!/bin/bash
#
# ===== NOTE ========
# This script launches a timing_code on a Xeon

# This is the executable. In my case the host directory is mounted on the Xeon Phi
# Your mileage may vary
EXECUTABLE=`pwd`/build_avx/tests/time_dslash_noqdp
EXECUTABLE_ARGS="-x 32 -y 32 -z 32 -t 64 -by 8 -bz 8 -pxy 1 -pxyz 0 -c 16 -sy 1 -sz 2 -prec f -i 300 -minct 2 -compress12"

# This is the location of the OMP Shared libraries on the Xeon Phi
# Your mileage may vary
export KMP_AFFINITY=compact
${EXECUTABLE} ${EXECUTABLE_ARGS}
