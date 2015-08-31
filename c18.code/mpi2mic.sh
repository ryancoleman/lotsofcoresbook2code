#!/bin/env bash

#copy bindaries in mic0:/tmp/
make miccopy

#copy libraries and pmi_proxy in mic0:/tmp/
make miclibcopy

export I_MPI_MIC=enable 

mpirun -n 1 -host mic0 -genv LD_LIBRARY_PATH /tmp \
-genv OMP_PLACES threads -genv OMP_PROC_BIND spread,close \
-genv OMP_NESTED true -genv OMP_NUM_THREADS 60,4 \
-genv KMP_HOT_TEAMS_MAX_LEVEL 2  -genv KMP_HOT_TEAMS_MODE 1 /tmp/dgemm

mpirun -n 1 -host mic0 -genv LD_LIBRARY_PATH /tmp \
-genv OMP_PLACES threads -genv OMP_PROC_BIND spread,close \
-genv OMP_NESTED true -genv OMP_NUM_THREADS 60,4 \
-genv KMP_HOT_TEAMS_MAX_LEVEL 2  -genv KMP_HOT_TEAMS_MODE 1 /tmp/fft3d -g 64 -i 10

mpirun -n 1 -host mic0 -genv LD_LIBRARY_PATH /tmp \
-genv OMP_PLACES threads -genv OMP_PROC_BIND close -genv OMP_NUM_THREADS 240 \
/tmp/diffusion mic nx=240 nf=10 out=
