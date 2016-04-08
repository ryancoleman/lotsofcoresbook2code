#!/bin/bash
#
# ===== NOTE ========
# This script launches a timing_code on a Xeon Phi, using an MPI launcher
# For this to work you must have built the Xeon Phi target in this directory using the ./build_mic.sh 
# script
#
# For this to be succesful, the MPI must be set up (mpivars.sh must be sourced)
# The executable has to be located on the XeonPhi card
# 
# Either the OpenMP shared libs must be on the LD_LIBRARY_PATH automatically or they need
# to be set explicitly

# This is the hostname of my target Xeon Phi
# Your mileage may vary
TARGET_MIC=qcd12m0409-mic0

# This is the executable. In my case the host directory is mounted on the Xeon Phi
# Your mileage may vary
EXECUTABLE=`pwd`/build_mic/tests/time_dslash
EXECUTABLE_ARGS=""

# This is the location of the OMP Shared libraries on the Xeon Phi
# Your mileage may vary
OMPLOCATION=/opt/intel/composerxe/lib/mic

export I_MPI_MIC=1
mpiexec.hydra -n 1 -host ${TARGET_MIC} -genv LD_LIBRARY_PATH=${OMPLOCATION}:$LD_LIBRARY_PATH -env KMP_AFFINITY=compact,granularity=thread -env OMP_NUM_THREADS=236 ${EXECUTABLE} ${EXECUTABLE_ARGS}

