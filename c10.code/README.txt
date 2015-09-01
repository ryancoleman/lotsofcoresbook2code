
# MODAL mini-app

## Overview
Part 1 contains 10 versions of the code, each one adding further optimisations to the code. The makefile builds an executable for each. Parts 1-3 have no offloading.
Part-2 explores different nested parallelism models and has an executable for each one.

## Requirements for build
* Intel icc and icpc versions 15+.
* Intel MKL
* GNU Scientific Library (GSL)
* iniparser (https://github.com/ndevilla/iniparser)

Build iniparser from source:
From base directory cd to resources/iniparser.
	Type `make`
	You may need to append [BASEDIR]/resources/iniparser to your
          LD_LIBRARY_PATH. 

## Building
Build part-1 with offloading with:
$ make part-1

Build part-1 without offloading with:
$ make part-1-no-offload

Build part-2 with offloading with:
$ make part-2

Build all the executables for part-1 and part-2 with offloading enabled with:
$ make all

Inside the Makefile is an option L1_CUTOFF. This shortens the problem size so that the runtime is short enough for benchmarking.
Typical value is 50.

## Running

To run executables from part-1 do:
$ source part-1/offload_vars.sh
$ mpirun -np 1 part-1/modal_V[n] [parameters.ini]

Where [n] is the version number which has values 1 to 10 and [parameters.ini] is one of the .ini files in the base directory.

To run executables from part-2, run one of the following scripts:
manual.sh   (uses manual nested parallelism)
omp.sh      (uses nested openmp regions)
omp_hot_teams.sh    (uses nested openmp regions with hot teams environment variables)
teams.sh    (uses openmp teams construct)
