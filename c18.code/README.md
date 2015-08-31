Test codes
================
- dgemm.cpp : MPI and Nested OpenMP runs using dgemm in MKL
- fft3d.cpp : MPI and Nested OpenMP runs using 3D FFT in MKL
- diffusion.c : Various OpenMP implementations

How to build
================

Makefile is provided to build the binaries for the host and those for KNC
in mic directory.

To compile with Intel(R) C/C++ Compiler for Linux,

* set the environments 

  $ source <path_to_compiler_bin>/iccvars.sh or iccvars.csh
  $ source <path_to_mpi_root>/bin64/mpivars.sh
 
 - Set ICC_ROOT and IMPI_ROOT in Makefile

  ICC_ROOT=/opt/intel/composer_xe_2015.2.164
  IMPI_ROOT=/opt/intel/impi/5.0.3.048/

* set HOSTFLAGS for the Host. Currently it is set to -xAVX 
  
  $ HOSTFLAGS   = -xAVX
 
* build binaries in host and mic directories
  $ make


To copy binaries and libraries on mic0:/tmp
  $ make miccopy
  $ make miclibcopy

The common Intel compiler options are set in ICFLAGS.

  ICFLAGS    =  -restrict -unroll  -O3 -qopenmp -g  -DINLINE=inline -DMUST_INLINE=inline -DHAVE_MKL

Note that the options are for Intel Composer 15.x and higher.

How to check 
=============

* Run on a host: make run_host
 - Execute host/dgemm, host/fft3d and host/diffusion

* Run on mic0 using ssh
 - make mic_dgemm
 - make mic_fft3d
 - make mic_diffusion

* Run on mic0 using mpi: use mpi2mic.sh script
 - source mpi2mic.sh

How to run on host
==================

1. dgemm.cpp : multiply_mkl calls cblas_dgemm 

  mpirun -np MPI_TASK ./dgemm [M N K]

  Defaults: M=N=K=256

  The number of iterations is set in the code based on the problem sizes.

2. fft3d.cpp 

  mpirun -np MPI_TASK ./fft3d -g 64 -i 100

  Options
  -g int : grid sizes for x, y and z
  -x int : grid size for x
  -y int : grid size for y
  -z int : grid size for z
  -i int : number of iterations
  
3. diffusion.c

  ./diffusion METHOD nx=256 count=1200 ouput=

  Choose METHOD among base, openmp, mic, nested, nested_lb, twoyz, task, crew
  Options
  nx=int    The x, y, and z dimensions, default=256
  count=int The number of iterations (per frame), default=1200
  nf=int    The number of frames, default=10
  out=xxx   Optional outpu file, default is none

How to run on a coprocessor
=============================

For comprehensive benchmarking, use these scripts on a corprocessor.
- dgemm.sh
- fft3d.sh
- diffusion.mic.sh

Three shell scripts for a series of tests are privded for each compute.  All
the tests can be excecuted *natively* (ssh mic0), as far as all MKL, MPI and
TBB libraries can be found in `LD_LIBRARY_PATH` and MPI runtime is in `PATH`.

Make sure LD_LIBRARY_PATH points to the directory with MPI, OpenMP and MKL
dynamic libraries and PATH to the directory where pmi_proxy can be found.

