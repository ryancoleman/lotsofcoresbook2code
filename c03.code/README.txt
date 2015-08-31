1. Type 'make' to compile WRF Goddard scheme for both CPU and Xeon Phi coprocessor

2. To execute CPU code, for both the original and optimized codes, run script 'run_CPU_code.sh'

3. To execute Xeon Phi codes run script 'run_xeon_phi_code.sh'
   -You might need to edit run_xeon_phi_code.sh to set the location of Xeon Phi's OpenMP library on your system (locate it using e.g. 'locate mic/libiomp5.so' command)


The are several options in the standalone driver, gsfcgce_driver.f90:

-Default data directory is controlled by this statement  
  call chdir("data") ! data directory for both input and output files

-You can set the number of execution runs using this parameter
  integer, parameter :: no_runs = 10

-The following parameter controls weather or not the outputs from a exection are written to output files
  logical, parameter :: write_output = .true.

