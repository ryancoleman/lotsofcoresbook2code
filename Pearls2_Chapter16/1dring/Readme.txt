The code for 1-d ring example, described in "An Introduction to MPI-3 Shared
Memory Programming" 


* How to build and run:
This code was tested with Intel® MPI library 5.0.2 and Intel® C++ Compiler v15.0.1 on Red Hat* Enterprise Linux 6.5 OS.
Note: compiler and MPI environments should be set, like follows (probably your pathes will be different)
. /opt/compilers/intel/composer_xe_2015.2.164/bin/compilervars.sh intel64
. /opt/intel/impi/5.0.2.044/bin64/mpivars.sh 
 


Build:
- make       # builds for both host and mic 
- make host  # builds for host
- make mic   # builds for mic


Run:
Note: runs  with coprocessor will work if you have file system shared with corpocessor and the binary is available for coprocessor. 
Otherwise see howto instructions in the end of this file 

To run the program with 2 processes residing on the host only the following command may be used:
- make run-host

To run the program with 2 processes residing on the coprocessor only the following command may be used:
- make run-mic

To run the program with 4 processes (2 of them on the host and 2 are on coprocessor) the following command may be used:
- make run-symm

To run the program with more processes and/or on different nodes "mpirun" command (which is typically a part of MPI distribution) should be used. 
An example command line which runs 16 processes (4 processes per node)  is:

mpirun -n 16 -ppn 4 -f hostfile.txt ./mpi3shm_1Dring.v2.1
where hostfile contains host names (1 per line), like:
hostname01
hostname02
hostname03
hostname04


If there is no shared file system with coprocessor the following steps should
be taken to run the programm on coprocessor only:
- copy the binary with .mic prefix to the coprocessor (via scp for instance)
- perform the following command:

I_MPI_MIC=1 mpirun -n 16 -host hostname-mic0 ./mpi3shm_1Dring.v2.1.mic
 where "hostname" is the name of your machine

To run the programm in symmetric mode additional steps are required:  path to
the binary on coprocessor should be specified as well as  .mic extension. So,
the command will look like follows:

I_MPI_MIC=1 mpirun -n 16 -ppn 8 -hosts hostname,hostname-mic0 -genv I_MPI_MIC_PREFIX path_to_binary -genv I_MPI_MIC_POSTFIX .mic  ./mpi3shm_1Dring.v2.1
where "hostname" is the name of your machine and "path_to_binary" is the folder containing the binary on coprocessor
   
 
  
 
