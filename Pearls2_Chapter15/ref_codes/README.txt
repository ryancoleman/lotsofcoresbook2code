There are 3 different application programs in this directory: Matrix
Multiply, Cholesky Factorization, and LU Factorization. These are 
contained in the sub-directories matMult, cholesky, and lu, respectively. 

For Cholesky and LU programs, two versions are available, one for host 
(CPU) only code (named "tiled_host"), and other named "tiled_hstreams" 
which corresponds to an offload code utilizing the Xeon Phi card as well 
as the host CPU through the hStreams framework.

There are separate README files present in each subdirectory for the 
three applications which describe in detail the program as well as the
build and run instructions. 

There is a top-level Makefile present which builds all the programs.
Also, there are individual makefiles in the three subdirectories for
matMult, cholesky, and lu, which build individual applications.
