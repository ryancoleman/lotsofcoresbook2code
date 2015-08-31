#!/bin/bash 

mpiifort -openmp -g -traceback -O3 -o stest_hybrid_part_coarse.x \
         stest_hybrid_part_coarse.f90  \
         shado_xchange_EW_part_coarse.f90  \
         shado_xchange_NS_part_coarse.f90  \
         shado_xchange_TB_part_coarse.f90

# These values are just placeholders, change as required:
export OMP_NUM_THREADS=4
export I_MPI_PIN_DOMAIN=omp
