The updated MPPTEST code, which includes MPI-3 SHM halo patterns, described in "An Introduction to MPI-3 Shared Memory Programming"  chapter


* How to build and run:
This code was tested with the following SW and HW configuration:

- dual Intel® Xeon® E5-2697 processors with one Intel® Xeon Phi™ 7120P-C0 coprocessor
  and one Mellanox* Connectx-3 InfiniBand* adapter connected to the same socket. 
- Red Hat* Enterprise Linux 6.5 OS, 
- Intel MPSS 3.3.30726 
- OFED 1.5.4.1 
- Intel® MPI Library v5.0.2 
- Intel® C++ Compiler v15.0.1


Build:

MPI distribution is required to build mpptest 

The typical configuration string is below. 
(you will probably have to specify another desctination path and path to mpi distribution relevant for your system). 
Note that MPICC needs to refer to mpiicc, which uses Intel compiler.    

export MPICC="mpiicc"

./configure --with-mpi=/opt/intel/impi.5.0.2/intel64/ CFLAGS="-DHAVE_MPI_PASSIVERMA -DHAVE_MPI_GET -DHAVE_MPI_PUT"

make (and "make install" if there is a need)

To build the program for many core arcitecture, such as Intel® Xeon Phi™ coprocessor the following variable should be set prior to configure:

export MPICC="mpiicc -mmic"

Then the same configure string can be used. Note that MPSS is required for building for MIC. 
   
./configure --with-mpi=/opt/intel/impi.5.0.2/intel64/ CFLAGS="-DHAVE_MPI_PASSIVERMA -DHAVE_MPI_GET -DHAVE_MPI_PUT"

make EXEEXT=.mic (and "make install" if there is a need)

This will produce "mpptest.mic" binary compiled for Intel® Xeon Phi™ coprocessor.

Note that "make clean" could be needed if build for other platform has been done previously.

Run:

run_mpi3shm_mpptest.sh file may be used as a run script. It describes all
needed actions in details. Please note that you need to create hostfile
specific to your system (see example in run_mpi3shm_mpptest.sh)      
Below are some hints:


- there should be shared file system with coprocessor to run the test in native or symmetric mode 
  (i. e. when all or some of the ranks residing on coprocessor). Otherwise the mpptest.mic binary 
  should be copied to coprocessor (by using scp for instance). How to setup NFS is described in 
  ch 6.4 "Setup Mounted File System" of MPSS CLuster Setup Guide, available here: 
  https://software.intel.com/en-us/articles/intel-manycore-platform-software-stack-mpss#downloads    

- for symmetric runs, where some ranks are on the host and others are on
  coprocessor the following environment variables may need to be used:
        I_MPI_MIC_POSTFIX=.mic  - specifies MIC binary suffix   
        I_MPI_MIC_PREFIX=path_to_mic_binary - specifies the path to binary on coprocessor (typically needed when 
        there is no shared file system and binaries are located in different folders on host and coprocessor). 



 
