1. Introduction

The source code in this package may be built and run on either
a standard Xeon host or a Xeon Phi coprocessor, and run using 
different permutations of OpenMP threads and MPI processes. 

The build and run scripts provide use the Intel Fortran compiler and Intel
MPI.  However, they can be easily modified to use other compilers or MPIs.

Corresponding to the terminology used in the accompanying book chapter,
there are 3 sub-directories here:
 fine-grain 
 partial-coarse-grain
 full-coarse-grain  

2. Building

The code in each directory may be built by typing "make" in each one, 
or by running one of the "build*.sh" scripts contained therein.

3. Running

See the run*.sh scripts in each sub-directory, which may be used as 
templates to run each case.  Sample pbs scripts are also provided. 
All such files will need to be edited to suit your own local installation.
In particular, please edit the LD_LIBRARY_PATH and MIC_LD_LIBRARY_PATH 
environment variables to list all directories containing libraries needed
at run-time (esp. compiler and MPI libraries).

By default, each executable reads its input data from "standard input",
but it is nearly always more practical to put the input data into a separate
file, and "re-direct" input from that file at run-time.  In fact the original
"prompts" for the input values have been commented out in the source, so 
unless provided as described here, the program will appear to "hang" at 
run-time!! 

First, then, prepare the "input" file which contains the "control" parameters,
and is read by the executable at the start of each run.  A sample input file
is provided in this directory (input_sample_6GB.txt), and may be modified as
desired. This file contains:

 40             # Number of smoothing iterations to be done.
 800, 800, 1200 # Number of points in each dimension of 3-d grid
 1              # Interval between smoothing iterations for screen output
 0.05           # Weight to be given to each neighbouring point (.001 - .1)
 100.           # Max. range of initial data (1. - 100000.)
 0              # Deterministic initialization if zero, random init otherwise
 2, 2, 6        # nprocx, nprocy, nprocz domain decomposition

The first parameter above (40) is the number of iterations to be performed. 
Increase this value for longer runs, or reduce it for shorter runs.

The parameters on the second line above (800, 800, 1200) define the number
of grid-points to be used in each dimemsion of the 3-d grid, and determine
the overall problem size.  These particular values provided correspond to 
total memory consumption of approx. 6GB.  Increase these numbers to create
a larger problem size, or reduce them for a smaller problem.  However, 
make sure that the "domain decomposition" parameters on the last line 
(2,2,6 in the file provided) divide evenly into the corresponding grid-size.

The parameter on the third line (1) determines the frequence for screen output
(i.e., 1 means every single iteration).  Usually leave this unchanged.

The parameter on the fourth and fifth lines (0.05, and 100., resp.) should
also be left unchanged, at least for benchmarking purposes. 

The parameter on the 6th line (0) selects a deterministic initialization
when zero (0), or a random initialization otherwise.  Deterministic 
initialization also means deterministic outut, so leave this at 0 to 
check that the output remain the same, e.g., as the number of threads 
or MPI processes are varied from run to run.

The values on the 7th (last) line above control the MPI "domain decomposition", 
and in particular the number of MPI processes to be used in the x, y and
z directions, respectively, of the regular 3-d "box" defined by the 
grid-points on line 2. Take care that nprocx divides evenly into the 
number of x-dimension grid-points; ditto for the y and z dimensions.  E.g., 
the (2,2,6) decomposition uses a total of 2x2x6=24 MPI processes. 

Output:  Note that output is written both to "standard output" (which can be 
re-directed as desired), and to a "/tmp/stresstest.dat" file (intended to 
provide a more permanent output record in case std. out. is lost).  Please 
check /tmp for this file.  Also beware that it can be over-written by 
subsequent runs. 

3.1 To run on a standard xeon host node

Assuming that an input file has been prepared (in the same directory as
the executable), then a job can be launched on a standard host using, e.g.:

export KMP_STACKSIZE=1g
export KMP_AFFINITY="scatter,granularity=fine"
export OMP_NUM_THREADS=2
export I_MPI_PIN_DOMAIN=omp

(time mpiexec.hydra -n 24 ./stest_hybrid_fine.x < input_sample_6GB.txt ) >& p24_t2_6GB.log

This runs a job using 24 MPI processes and 2 OpenMP threads per process.
Note the input redirection from the file input_sample_6GB.txt.  

3.2 To run on a Xeon Phi coprocessor

The basic settings and commands are similar to those above, with some extra 
environment variables for the Xeon Phi:

# These values are just placeholders, change as required:
export I_MPI_MIC=enable
export I_MPI_MIC_POSTFIX=.MIC
export I_MPI_PIN_DOMAIN=omp

export MIC_ENV_PREFIX=MIC_
export KMP_STACKSIZE=200MB
export MIC_KMP_STACKSIZE=50MB
export MIC_KMP_MONITOR_STACKSIZE=12MB
export KMP_AFFINITY="scatter,granularity=fine"
export MIC_KMP_AFFINITY="scatter,granularity=fine"

export OMP_NUM_THREADS=100
export MIC_OMP_NUM_THREADS=100

(time mpiexec.hydra -n 2 -hosts mic0 ./stest_hybrid_fine.x < input_sample_6GB_p2_t100.txt ) >& p2_t100_6GB_fine.log

Alternatively, "-machinefile ./hostfile" could be used instead of "-hosts mic0"
above, especially if you want to run on more than one MIC, or more than one 
node.  Here, "hostfile" contains the names of the nodes to be used, with a 
-mic suffix, e.g., perhaps just the single line:
phinally-mic0

There are many variations on how such a job can be launched on Phi coprocessors.
The key points are to make the thread-count and MPI decomposition consistent 
with the parameter settings in the input file.

Further environment variables may need to be defined, e.g., LD_LIBRARY_PATH 
and MIC_LID_LIBRARY_PATH, so that all run-time libraries needed by the 
executables can be found.  In particular, if MIC_ENV_PREFIX=MIC_ is used,
then the "mic" compiler libraries should be appended to LD_LIBRARY_PATH 
(as is done in the run*sh templates provided).

Enda O'Brien
enda.obrien@ichec.ie

