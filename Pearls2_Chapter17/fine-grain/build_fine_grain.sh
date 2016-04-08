#!/bin/bash 

mpiifort -openmp -g -traceback -O3 -o stest_hybrid_fine.x \
         stest_hybrid_fine.f  \
         shado_xchange_EW_fine.f  \
         shado_xchange_NS_fine.f  \
         shado_xchange_TB_fine.f

# These values are just placeholders, change as required:
export OMP_NUM_THREADS=4
export I_MPI_PIN_DOMAIN=omp
