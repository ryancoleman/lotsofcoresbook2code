# FILE LICENSE TAG: SAMPLE 

README for LU Factorization Reference Code (2nd February, 2015) 

This file is for use for the LU on Linux only.

For information on windows, please see the README.txt file in ../windows/lu.

This file contains the following sections:

Contents
* INTRODUCTION
* HOW TO BUILD THE LU REFERENCE CODE
* HOW TO RUN THE LU REFERENCE CODE

INTRODUCTION

LU decomposition (or LU factorization), where 'LU' stands for 'Lower Upper',
is a decomposition of a matrix into the product of a lower triangular matrix 
and a upper triangular matrix. It is useful for efficient numerical solutions
and Monte Carlo simulations. The LU decomposition is just a matrix form of
Gaussian elimination [Ref. 1].

Let [A] be a square matrix having M rows by M columns. (The algorithm can also
be applied to rectangular matrices, but this is rarely done in practice, and 
the supplied reference code only works for square matrices). The idea is to 
transform [A] into an M by M upper triangular matrix [U] by introducing zeros 
below the diagonal, first in column 1, then in column 2, and so on. This is done
by subtracting multiples of each row from subsequent rows. This "elimination"
process is equivalent to multiplying [A] by a sequence of lower-triangular 
matrices [L_k] on the left:
    [L_m-1]...[L_2][L_1][A] = [U]
and setting [L] = [L_1]^(-1)[L_2]^(-1)...[L_m-1]^(-1), gives
    [A] = [L][U]
Thus we obtain LU factorization of [A], where [U] is upper-triangular and [L]
is unit lower-triangular (which means all of its diagonal entries are equal
to 1) [Ref. 1]. 

The computational cost of a LU decomposition is (2/3)*M^3 floating point
operations. All the three matrices [A], [L], and [U] are not needed to be stored;
to minimize memory use, [L] and [U] are overwritten on input matrix [A]. 
Thus the memory footprint is M^2*sizeof(double) bytes (for double-precision real-
valued matrix). An algorithm which is applied on elements of matrix [A] can 
only proceed one column at a time as there is a data dependence which limits 
the overall concurrency. The reference code implements a right-looking LU
decomposition algorithm and therefore as we proceed down the columns, the 
number of entries updated reduces [Ref. 1].

In this program, we do a tiled implementation of LU decomposition without 
pivoting (for double-precision real-valued matrices) using fast level-3 BLAS
components (using the  Intel(R) MKL library).  The tiled-algorithm is similar 
to the one for Cholesky factorization as given in Ref. 2. We use the hStreams 
framework to partition the card into 4 or 6 physical partitions and schedule 
BLAS3 operations (DGETRF, DTRSM, DGEMM) into compute-streams associated with 
these partitions in an as concurrent manner as possible. Note that for DGETRF, 
we use a custom hand-written kernel instead of the one from MKL. This is because 
we implement the LU decomposition without pivoting, whereas the function available
in MKL does pivoting. The matrix is divided into T x T tiles each containing 
(M/T * M/T elements). The outer-loop is iterated over 1:T. For the k-th outer-
loop iteration, we have 1 DGETRF (of the k x k tile), 2*(T-k) DTRSMs, and (T-k)*(T-k) 
DGEMMs. The order of operations (and dependence) of BLAS3 operations for a given k 
is: DGETRF -> DTRSM -> DGEMM (i.e. DGEMM depends on DTRSM, which in turn depends 
on DGETRF). Thus the max concurrency for a given k is : no. of DGEMMs = (T-k)*(T-k). 
The concurrency is highest for k = 1, and is = (T-1)*(T-1). If T = 6, for 
example, then the max concurrency is 25. It's best if this max concurrency is 
divisible by number of partitions on the card. As a rule of thumb, 
we have observed that no. of partitions on card = T - 1 gives the best results, 
and in particular, we have observed that T = 6, and no. of partitions = 5, 
give best results for most cases.

There are two versions of the tiled-LU program. One called tiled_host performs 
the tiled-LU on the host only (without using hStreams or the card). The 
performance of this program is compared against the MKL DGETRF without 
automatic offload. The other version is called tiled_hstreams. This
uses both the host and the card, by using the hStreams library. The performance
of the tiled_hstreams version is compared against the MKL DGETRF with automatic
offload. For the tiled_hstreams LU decomposition, data transfer to/from the 
card is interleaved with compute on the card using the concurrently running streams. 

To understand the data/compute dependencies, it's recommended that the user
first becomes familiar with the algorithm in the tiled_host example and then
it will be easier to understand the various synchronizations (in the form of 
_event_wait) required in the tiled_hstreams example.

References:
1) Trefethen, Lloyd N. and Bau III, David, Numerical Linear Algebra, SIAM (1997).
2) Jeannot, Emmanuel. Performance Analysis and Optimization of the Tiled 
Cholesky Factorization on NUMA Machines. PAAP 2012-IEEE International Symposium
on Parallel Architectures, Algorithms and Programming. 2012.

In the description above, square brackets indicate a matrix
(e.g. [A]).  The symbol ^ indicates exponentiation by a power. The
symbol * indicates multiplication.


1. To run the tiled_host example
You would run the lu_tiled_host executable as follows:

./lu_tiled_host -m 4800 -t 6 -l row -i 5

the description of command line inputs is:
-m : matrix size per side
-t : no. of tiles the matrix is broken into per side
-l : row (or ROW) if ROW_MAJOR layout of matrix storage is used, otherwise COL_MAJOR
-i : no. of iterations to be performed

NOTE: the matrix size must be divisible by number of tiles

one can also use the runit.sh script provided, which also sets the KMP_AFFINITY. 

2. To run the tiled_hstreams example
You would run the lu_tiled_hstreams executable as follows:

./lu_tiled_hstreams -m 4800 -t 6 -s 5 -l row -i 5

the description of command line inputs is:
-m : matrix size per side
-t : no. of tiles the matrix is broken into per side
-s : no. of logical streams (= no. of physical partitions on the card)
-l : row (or ROW) if ROW_MAJOR layout of matrix storage is used, otherwise COL_MAJOR
-i : no. of iterations to be performed

NOTE: the matrix size must be divisible by number of tiles

one can also use the runit.sh script provided, which also sets the KMP_AFFINITY on host and card. 


HOW TO BUILD THE LU REFERENCE CODE

1. Install MPSS 3.5
2. Install the Intel Composer XE compiler
3. Copy the reference code to an empty temporary directory:
   $ cd
   $ rm -fr temp_ref_code
   $ mkdir temp_ref_code
   $ cd temp_ref_code
   $ cp -r /usr/share/doc/hStreams/ref_code .
4. Change directory to the ref_code/lu dir
   $ cd ref_code/lu
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
icc       -openmp  -I../../common -mkl -DNDEBUG -O3 -Wl,--enable-new-dtags  -c lu_tile.cpp -o lu_tile.host.o
icc       -openmp  -I../../common -mkl -DNDEBUG -O3 -Wl,--enable-new-dtags  -o ../../../bin/host/lu_tiled_host lu_tile.host.o ../../common/dtime.host.o ../../common/matrices_generator.host.o

7. Now build the tiled_hstreams version. Change directory to the ../tiled_hstreams dir
   $ cd ../tiled_hstreams
8. Type make:
   $ make

The make should build as follows:
$ make
icc       -openmp -I/usr/include/intel-coi -I../../common -I../../../include -I../../../src/include -mkl -DNDEBUG -O3  -Wl,--enable-new-dtags  -lcoi_host -c lu_tile_hstreams.cpp -o lu_tile_hstreams.host.o
icc       -openmp -I/usr/include/intel-coi -I../../common -I../../../include -I../../../src/include -mkl -DNDEBUG -O3  -Wl,--enable-new-dtags  -lcoi_host -c hStreams_custom.cpp -o hStreams_custom.host.o
icc       -openmp -I/usr/include/intel-coi -I../../common -I../../../include -I../../../src/include -mkl -DNDEBUG -O3  -Wl,--enable-new-dtags  -lcoi_host -o ../../../bin/host/lu_tiled_hstreams lu_tile_hstreams.host.o ../../common/dtime.host.o hStreams_custom.host.o ../../common/matrices_generator.host.o -L../../../bin/host  -lhstreams_source
icc -mmic -openmp -I/usr/include/intel-coi -I../../common -I../../../include -I../../../src/include -mkl -DNDEBUG -O3  -Wl,--enable-new-dtags  -lcoi_device -rdynamic -fPIC -shared -c hStreams_custom_kernels_sink.cpp -o hStreams_custom_kernels_sink.dev.o
icc -mmic -openmp -I/usr/include/intel-coi -I../../common -I../../../include -I../../../src/include -mkl -DNDEBUG -O3  -Wl,--enable-new-dtags  -lcoi_device -rdynamic -fPIC -shared -o ../../../bin/dev/lu_sink_1.so hStreams_custom_kernels_sink.dev.o -L../../../bin/dev -lhstreams_sink -Wl,-soname,lu_sink_1.so -static-intel
icc: warning #10342: -liomp5 linked in dynamically, static library not available for Intel(R) MIC Architecture
icc: warning #10342: -liomp5 linked in dynamically, static library not available for Intel(R) MIC Architecture
icc: warning #10237: -lcilkrts linked in dynamically, static library not available


HOW TO RUN THE LU REFERENCE CODE

Run the tiled_host version
1. Change directory to the tiled_host dir, and run the runit.sh script:

   $ cd ../tiled_host
   $ ./runit.sh 

This is how a typical output looks:

matrix is in Row major format
mat_size = 4800, num_tiles = 6, niter = 5


iter 0
time for Tiled LU = 521.655083 msec
time for MKL LU (dgetrf) (host, NO AO) = 383.312941 msec
Tiled LU successful

iter 1
time for Tiled LU = 486.261845 msec
time for MKL LU (dgetrf) (host, NO AO) = 350.229979 msec
Tiled LU successful

iter 2
time for Tiled LU = 479.424953 msec
time for MKL LU (dgetrf) (host, NO AO) = 377.222061 msec
Tiled LU successful

iter 3
time for Tiled LU = 480.629921 msec
time for MKL LU (dgetrf) (host, NO AO) = 382.603168 msec
Tiled LU successful

iter 4
time for Tiled LU = 479.327917 msec
time for MKL LU (dgetrf) (host, NO AO) = 366.101980 msec
Tiled LU successful

Tiled LU, for 4 iterations (ignoring first),
mean Time = 481.41 msec, stdDev Time = 3.29 msec,
Mean Gflops (using mean Time) = 153.15

MKL LU, For 4 iterations (ignoring first),
mean Time = 369.04 msec, stdDev Time = 14.30 msec,
Mean Gflops (using mean Time) = 199.78


Run the tiled_hstreams version
1. Change directory to the tiled_hstreams dir, and run the runit.sh script:

   $ cd ../tiled_hstreams
   $ ./runit.sh 

This is how a typical output looks:

~/hstreams/bin/host ~/hstreams/ref_code/lu/tiled_hstreams
~/hstreams/bin/dev ~/hstreams/bin/host ~/hstreams/ref_code/lu/tiled_hstreams
~/hstreams/bin/host ~/hstreams/ref_code/lu/tiled_hstreams
~/hstreams/ref_code/lu/tiled_hstreams
matrix is in Row major format
no. of streams (partitions) = 5, mat_size = 4800, num_tiles = 6, niter = 5


Iteration = 0, Tbegin
time for Tiled hstreams LU for iteration 0 = 818.50 msec
time for MKL LU (dgetrf) (AO) for iteration 0 = 276.64 msec
Tiled LU successful

Iteration = 1, Tbegin
time for Tiled hstreams LU for iteration 1 = 336.69 msec
time for MKL LU (dgetrf) (AO) for iteration 1 = 334.06 msec
Tiled LU successful

Iteration = 2, Tbegin
time for Tiled hstreams LU for iteration 2 = 388.07 msec
time for MKL LU (dgetrf) (AO) for iteration 2 = 382.38 msec
Tiled LU successful

Iteration = 3, Tbegin
time for Tiled hstreams LU for iteration 3 = 346.46 msec
time for MKL LU (dgetrf) (AO) for iteration 3 = 360.69 msec
Tiled LU successful

Iteration = 4, Tbegin
time for Tiled hstreams LU for iteration 4 = 471.69 msec
time for MKL LU (dgetrf) (AO) for iteration 4 = 377.19 msec
Tiled LU successful

Matrix size = 4800
Tiled hStreams LU: for 4 iterations (ignoring first),
mean Time = 385.73 msec, stdDev Time = 61.49 msec,
Mean Gflops (using mean Time) = 191.14

MKL AO LU (dgetrf): For 4 iterations (ignoring first),
mean Time = 363.58 msec, stdDev Time = 21.74 msec,
Mean Gflops (using mean Time) = 202.78

