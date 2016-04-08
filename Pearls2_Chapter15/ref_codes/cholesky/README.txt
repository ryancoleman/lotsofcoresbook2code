# FILE LICENSE TAG: SAMPLE 

README for Cholesky Factorization Reference Code (16th January, 2015) 

This file is for use for the cholesky on Linux only.

For information on windows, please see the README.txt file in ../windows/cholesky.

This file contains the following sections:

Contents
* INTRODUCTION
* HOW TO BUILD THE CHOLESKY REFERENCE CODE
* HOW TO RUN THE CHOLESKY REFERENCE CODE

INTRODUCTION

Cholesky decomposition is a decomposition of Hermitian, postive-definite matrix
into the product of a lower triangular matrix and its conjugate transpose, 
useful for efficient numerical solutions and Monte Carlo simulations. When it 
is applicable, the Cholesky decomposition is roughly twice as efficient as LU 
decomposition for solving systems of linear equations [Ref. 1].

If we define a square matrix [A] having M rows by M columns, then the Cholesky
decomposition is of the form:
    [A] = [L] * [L^t] 
where, [L] is a lower triangular matrix, and [L^t] is its transpose.
Every Hermitian positive-definite matrix (and thus also every real-valued 
symmetric postive-definite matrix) has a unique Cholesky decomposition.

The computational cost of a Cholesky decomposition is (1/3)*M^3 floating point
operations. Typically only half of the matrix A is stored (lower triangular),
and this matrix is overwritten by the cholesky factor matrix. Thus the memory 
footprint is (1/2)*M^2*sizeof(double) bytes (for double-precision real-valued 
matrix)  (In the reference code we do save the complete matrix).  An algorithm 
which is applied on elements of matrix [A] can only proceed one column at a 
time as there is a data dependence which limits the overall concurrency. 
Moreover, as one proceeds down the columns, the number of entries updated 
reduces as we have implemented the right-looking Cholesky decomposition algorithm
[Ref. 2].

In this program, we do a tiled implementation of Cholesky decomposition (for 
double-precision real-valued matrices) [Ref. 3] using fast level-3 BLAS
components. We use the hStreams framework to partition the card into 4 or 6 
physical partitions and schedule BLAS3 operations (DPOTRF, DTRSM, DSYRK, DGEMM)
into compute-streams associated with these partitions in an as concurrent 
manner as possible. The matrix is divided into T x T tiles each containing 
(M/T * M/T elements). The outer-loop is iterated over 1:T. For the k-th outer-
loop iteration, we have 1 DPOTRF (of the k x k tile), T-k DTRSMs, T-k DSYRKs,
and (1/2)*(T-k-1)*(T-k) DGEMMs. The order of operations (and dependence) of BLAS3
operations for a given k is: DPOTRF -> DTRSM -> DSYRK + DGEMM (i.e. DSYRK and 
DGEMM depend on DTRSM, which in turn depends on DPOTRF). Thus the max 
concurrency for a given k is : no. of DGEMMs + no. of DSYRKs = T-k + 
(1/2)*(T-k-1)*(T-k). The concurrency is highest for k = 1, and is = T-1 + 
(1/2)*(T-2)*(T-1). If T = 6, for example, then the max concurrency is 5 + 10 = 15.
It's best if this max concurrency is divisible by number of partitions on the
card. As a rule of thumb, we have observed that no. of partitions on card = 
T - 1 gives the best results, and in particular, we have observed that T = 6, 
and no. of partitions = 5, give best results for most cases.

There are two versions of the tiled-Cholesky program. One called tiled_host 
performs the tiled-Cholesky on the host only (without using hStreams or the 
card). The performance of this program is compared against the MKL DPOTRF 
without automatic offload. The other version is called tiled_hstreams. This
uses both the host and the card, by using the hStreams library. The performance
of the tiled_hstreams version is compared against the MKL DPOTRF with automatic
offload. For the tiled_hstreams Cholesky decomposition, data transfer 
to/from the card is interleaved with compute on the card using the concurrently
running streams. 

To understand the data/compute dependencies, it's recommended that the user
first becomes familiar with the algorithm in the tiled_host example and then
it will be easier to understand the various synchronizations (in the form of 
_event_wait) required in the tiled_hstreams example.

References:
1) Wikipedia. http://en.wikipedia.org/wiki/Cholesky_decomposition
2) Trefethen, Lloyd N. and Bau III, David, Numerical Linear Algebra, SIAM (1997).
3) Jeannot, Emmanuel. Performance Analysis and Optimization of the Tiled 
Cholesky Factorization on NUMA Machines. PAAP 2012-IEEE International Symposium
on Parallel Architectures, Algorithms and Programming. 2012.

In the description above, square brackets indicate a matrix
(e.g. [A]).  The symbol ^ indicates exponentiation by a power. The
symbol * indicates multiplication.


1. To run the tiled_host example
You would run the cholesky_tiled_host executable as follows:

./cholesky_tiled_host -m 4800 -t 6 -l col -i 5

the description of command line inputs is:
-m : matrix size per side
-t : no. of tiles the matrix is broken into per side
-l : row (or ROW) if ROW_MAJOR layout of matrix storage is used, otherwise COL_MAJOR
-i : no. of iterations to be performed

NOTE: the matrix size MUST be divisible by number of tiles

one can also use the runit.sh script provided, which also sets the KMP_AFFINITY. 

2. To run the tiled_hstreams example
You would run the cholesky_tiled_hstreams executable as follows:

./cholesky_tiled_hstreams -m 4800 -t 6 -s 5 -l col -i 5

the description of command line inputs is:
-m : matrix size per side
-t : no. of tiles the matrix is broken into per side
-s : no. of logical streams (= no. of physical partitions on the card)
-l : row (or ROW) if ROW_MAJOR layout of matrix storage is used, otherwise COL_MAJOR
-i : no. of iterations to be performed

NOTE: the matrix size MUST be divisible by number of tiles

one can also use the runit.sh script provided, which also sets the KMP_AFFINITY on host and card. 


HOW TO BUILD THE CHOLESKY REFERENCE CODE

1. Install MPSS 3.5
2. Install the Intel Composer XE compiler
3. Copy the reference code to an empty temporary directory:
   $ cd
   $ rm -fr temp_ref_code
   $ mkdir temp_ref_code
   $ cd temp_ref_code
   $ cp -r /usr/share/doc/hStreams/ref_code .
4. Change directory to the ref_code/cholesky dir
   $ cd ref_code/cholesky
5. Set the environment variables for the Intel Composer XE compiler:  

For example:

. /opt/mpss_toolchains/composer/composer_xe_2013/bin/compilervars.sh intel64
or
. /opt/intel/composerxe/bin/compilervars.sh intel64

(Your mileage may vary.  For example you probably will not have the Intel Composer
 XE compiler installed in /opt/mpss_toolchains).

5. First build the tiled_host version. Change directory to the tiled_host dir
   $ cd tiled_host
6. Type make:
   $ make

The make should build as follows:
$ make
icc       -openmp -I../../common -mkl -DNDEBUG -O3 -Wl,--enable-new-dtags  -c cholesky_tile.cpp -o cholesky_tile.host.o
icc       -openmp -I../../common -mkl -DNDEBUG -O3 -Wl,--enable-new-dtags  -o ../../../bin/host/cholesky_tiled_host cholesky_tile.host.o ../../common/dtime.host.o ../../common/matrices_generator.host.o

7. Now build the tiled_hstreams version. Change directory to the ../tiled_hstreams dir
   $ cd ../tiled_hstreams
8. Type make:
   $ make

The make should build as follows:
$ make
icc       -openmp -I/usr/include/intel-coi -I../../common -I../../../include -I../../../src/include -mkl -DNDEBUG -O3 -Wl,--enable-new-dtags  -lcoi_host -c cholesky_tile_hstreams.cpp -o cholesky_tile_hstreams.host.o
icc       -openmp -I/usr/include/intel-coi -I../../common -I../../../include -I../../../src/include -mkl -DNDEBUG -O3 -Wl,--enable-new-dtags  -lcoi_host -c hStreams_custom.cpp -o hStreams_custom.host.o
icc       -openmp -I/usr/include/intel-coi -I../../common -I../../../include -I../../../src/include -mkl -DNDEBUG -O3 -Wl,--enable-new-dtags  -lcoi_host -o ../../../bin/host/cholesky_tiled_hstreams cholesky_tile_hstreams.host.o ../../common/dtime.host.o hStreams_custom.host.o ../../common/matrices_generator.host.o -L../../../bin/host  -lhstreams_source
icc -mmic -openmp -I/usr/include/intel-coi -I../../common -I../../../include -I../../../src/include -mkl -DNDEBUG -O3 -Wl,--enable-new-dtags  -lcoi_device -rdynamic -fPIC -shared -c hStreams_custom_kernels_sink.cpp -o hStreams_custom_kernels_sink.dev.o
icc -mmic -openmp -I/usr/include/intel-coi -I../../common -I../../../include -I../../../src/include -mkl -DNDEBUG -O3 -Wl,--enable-new-dtags  -lcoi_device -rdynamic -fPIC -shared -o ../../../bin/dev/cholesky_sink_1.so hStreams_custom_kernels_sink.dev.o -L../../../bin/dev -lhstreams_sink -Wl,-soname,cholesky_sink_1.so -static-intel
icc: warning #10342: -liomp5 linked in dynamically, static library not available for Intel(R) MIC Architecture
icc: warning #10342: -liomp5 linked in dynamically, static library not available for Intel(R) MIC Architecture
icc: warning #10237: -lcilkrts linked in dynamically, static library not available


HOW TO RUN THE Cholesky REFERENCE CODE


Run the tiled_host version
1. Change directory to the tiled_host dir, and run the runit.sh script:

   $ cd ../tiled_host
   $ ./runit.sh 

This is how a typical output looks:

matrix is in Row major format
mat_size = 4800, num_tiles = 6, niter = 5


iter 0
time for Tiled Cholesky = 217.59 msec
time for MKL Cholesky (host, NO AO) = 328.01 msec
Tiled Cholesky successful

iter 1
time for Tiled Cholesky = 185.29 msec
time for MKL Cholesky (host, NO AO) = 316.62 msec
Tiled Cholesky successful

iter 2
time for Tiled Cholesky = 184.49 msec
time for MKL Cholesky (host, NO AO) = 320.98 msec
Tiled Cholesky successful

iter 3
time for Tiled Cholesky = 189.43 msec
time for MKL Cholesky (host, NO AO) = 328.06 msec
Tiled Cholesky successful

iter 4
time for Tiled Cholesky = 188.79 msec
time for MKL Cholesky (host, NO AO) = 336.56 msec
Tiled Cholesky successful

Tiled Cholesky, for 4 iterations (ignoring first),
mean Time = 187.00 msec, stdDev Time = 2.47 msec,
Mean Gflops (using mean Time) = 197.13

MKL Cholesky (host, NO AO), for 4 iterations (ignoring first),
mean Time = 325.55 msec, stdDev Time = 8.72 msec,
Mean Gflops (using mean Time) = 113.23


Run the tiled_hstreams version
1. Change directory to the tiled_hstreams dir, and run the runit.sh script:

   $ cd ../tiled_hstreams
   $ ./runit.sh 

This is how a typical output looks:

~/hstreams/bin/host ~/hstreams/ref_code/cholesky/tiled_hstreams
~/hstreams/bin/dev ~/hstreams/bin/host ~/hstreams/ref_code/cholesky/tiled_hstreams
~/hstreams/bin/host ~/hstreams/ref_code/cholesky/tiled_hstreams
~/hstreams/ref_code/cholesky/tiled_hstreams
matrix is in Row major format
no. of streams (partitions) = 5, mat_size = 4800, num_tiles = 6, niter = 5


Iteration = 0
time for Tiled hstreams Cholesky for iteration 0 = 394.65 msec
time for MKL Cholesky (AO) for iteration 0 = 290.59 msec
Tiled Cholesky successfull

Iteration = 1
time for Tiled hstreams Cholesky for iteration 1 = 176.76 msec
time for MKL Cholesky (AO) for iteration 1 = 325.66 msec
Tiled Cholesky successfull

Iteration = 2
time for Tiled hstreams Cholesky for iteration 2 = 170.63 msec
time for MKL Cholesky (AO) for iteration 2 = 330.30 msec
Tiled Cholesky successfull

Iteration = 3
time for Tiled hstreams Cholesky for iteration 3 = 180.76 msec
time for MKL Cholesky (AO) for iteration 3 = 314.84 msec
Tiled Cholesky successfull

Iteration = 4
time for Tiled hstreams Cholesky for iteration 4 = 171.93 msec
time for MKL Cholesky (AO) for iteration 4 = 336.41 msec
Tiled Cholesky successfull

Matrix size = 4800
Tiled hStreams Cholesky: for 4 iterations (ignoring first),
mean Time = 175.02 msec, stdDev Time = 4.65 msec,
Mean Gflops (using mean Time) = 210.62

MKL AO Cholesky: for 4 iterations (ignoring first),
mean Time = 326.80 msec, stdDev Time = 9.11 msec,
Mean Gflops (using meanTime) = 112.80

