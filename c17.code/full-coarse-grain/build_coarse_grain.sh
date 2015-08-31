#!/bin/bash 

mpiifort -openmp -g -traceback -mt_mpi -O3 -o stest_hybrid_coarse.x \
         stest_hybrid_coarse.f90  \
         shado_xchange_EW_coarse.f90  \
         shado_xchange_NS_coarse.f90  \
         shado_xchange_TB_coarse.f90

# These values are just placeholders, change as required:
export OMP_NUM_THREADS=4
export I_MPI_PIN_DOMAIN=omp
