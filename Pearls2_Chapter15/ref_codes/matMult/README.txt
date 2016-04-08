# FILE LICENSE TAG: SAMPLE 

README for matMult Reference Code (August 25, 2014) 

This file is for use of the matMult on Linux only.

For information on windows, please see the README.txt file in ../windows/matMult.

This file contains the following sections:

Contents
* INTRODUCTION
* HOW TO BUILD THE MATMUL REFERENCE CODE
* HOW TO RUN THE MATMUL REFERENCE CODE
* HOW TO INTERPRET RESULTS OF THE MATMUL REFERENCE CODE

INTRODUCTION

Matrix multiplication is an expensive mathematical operation used in
many software applications.  If we define two input matrices [A] and
[B] with [A] having M rows by K columns of numbers and [B] having K
rows by N columns of numbers, then the matrix product [C] is of size M
rows by N columns where each element of [C] is obtained by taking the
inner product (or dot product) of the corresponding row of [A] with
the corresponding column of [B].  The total number of multiplications
is (M*N*K) and the approximate number of additions is the same for a
total flops of (2*M*N*K).

If you try to use a high speed compute device to do matrix
multiplication, you need to send M*K*sizeof(double) bytes for matrix A
and then send K*N*sizeof(double) bytes for matrix B.  You typically
have to allocate enough memory (M*N*sizeof(double)) to hold the output
matrix first.  In certain common situations you have to initialize
that C output matrix to zero first. Then you start the matrix
multiplication operation on the high speed compute device.  And
finally you need to bring the output data (matrix C) back from the
device.

So the total data transfer is quantified as (M*K+K*N+M*N) double
precision numbers where double precision numbers require 8 bytes each.

If all three matrices in the expression [C] = [A]*[B] are square, then
the matrix multiplication operation requires 2*N^3 (2*Ncubed) Floating
Point Operations (flops) on the 24*N^2 (24*Nsquared) bytes.

If R_to is the data rate from the host to the compute device and
R_from is the data rate from the device back to the host, then the
"total data transfer time" T_xfer for square matrices is

       T_xfer = ((16*N^2 /R_to) + (8*N^2/ R_from)
         (bytes over bytes per second is seconds) 

If G_host is the achievable flops per second on the host, then the
"time to compute on the host" is

       T_host = (2*N^3 /G_host)
         (flops divided by flops  per second is seconds)

If G_card is the achievable flops on the card form factor compute
device, then the time to compute on the card is

       T_card = (2*N^3 /G_card )
         (flops divided by flops  per second is seconds)

To decide to use a separate card to compute with, one needs to decide
if

       T_card + T_xfer < T_host

Expanding the terms and approximating R_to = R_from ~=~ R yields the
following:

       (2*N^3 /G_card) + (24*N^2/R) < (2*N^3/G_card) 

Simplifying we get an equation with 3 reciprocals:

       (1/G_card) + (12/(R*N)) < (1/G_host) 

as the criterion for when to use the high speed compute device to do
the matrix multiplication.

Now what if we were able to asynchronously transfer data back and
forth to a card while asynchronously executing computations.  Then we
could block the [A] & [B] matrices into a set of blocks, multiply the
blocks on the high-speed-compute device, and then pull the results
back to the host.

To see the advantage, let F_row be a row blocking factor (e.g. 2, 3,
4, 5, etc) and let F_col be the column blocking factor.  Then the
first block of [A] (we will call it [A_00]) would be (M/F_row) x
(K/F_col) in size.  Similarly, the first block of [B] (we will call it
[B_00]) will be (K/F_row)x(N/F_col) in size.  This makes the initial
data transfer requires the following amount of time for square matrix
size N and square block factor F:

       T_blkxfer = 24*N^2/(R*F^2)

So the initial data transfer time can go down with the square of the
blocking factor.  If the asynchronous compute device is partitioned
into F partitions, we will have the asynchronous compute device
completely busy in the time proportional to (1/F).  Similarly, there
will be some time at the end of the operation where the device is not
completely busy and the last matrix block results are being copied off
of the device.  But overall, the net throughput of the asynchronous
compute device can be optimized by doing block matrix multiplication
rather than full size matrix multiplication if you can successfully
overlap data transfer and computation with low enough overhead for the
blocking.

So that is why we have looked at the performance of block matrix
multiplication using h-streams. To experiment with h-streams on Xeon
Phi(TM), we have provided an initial implementation that carries out
operations as described above.

To run the matMult example after defining the necessary environment
variables, you would use the mat_mult executable as follows:

       $ mat_mult -b1000 -m 8000 -n4000 -k16000 -i5 

will do an (8k by 16k) (MxK) matrix times a (16k by 4k) (KxN) matrix
to yield an (8k by 4k) (MxN) output matrix using 1k by 1k submatrix
blocks using internally 4 h-streams on the Xeon Phi.  The i5 option
tell the executable to do 5 iterations and leave off the first
iteration owing to initialization overhead.

In the description above, square brackets indicate a matrix
(e.g. [A]).  The symbol ^ indicates exponentiation by a power. The
symbol * indicates multiplication.


HOW TO BUILD THE MATMUL REFERENCE CODE

1. Install MPSS 3.5
2. Install the Intel Composer XE compiler
3. Copy the reference code to an empty temporary directory:
   $ cd
   $ rm -fr temp_ref_code
   $ mkdir temp_ref_code
   $ cd temp_ref_code
   $ cp -r /usr/share/doc/hStreams/ref_code .
4. Change directory to the ref_code/matMult dir
   $ cd ref_code/matMult
5. Set the environment variables for the Intel Composer XE compiler:  

For example:

. /opt/mpss_toolchains/composer/composer_xe_2013/bin/compilervars.sh intel64
or
. /opt/intel/composerxe/bin/compilervars.sh intel64

(Your mileage may vary.  For example you probably will not have the Intel Composer
 XE compiler installed in /opt/mpss_toolchains).

5. Type make:
   $ make

The make should build as follows:

$ make
icc       -openmp -I/usr/include/intel-coi -I../common -I../../include -I../../src/include -mkl -DNDEBUG -O3 -Wl,--enable-new-dtags -DFP_DATA_TYPE=double -lcoi_host -c matMultapi.cpp -o matMultapi.host.o
icc       -openmp -I/usr/include/intel-coi -I../common -I../../include -I../../src/include -mkl -DNDEBUG -O3 -Wl,--enable-new-dtags -DFP_DATA_TYPE=double -lcoi_host -o ../../bin/host/mat_mult matMultapi.host.o ../common/dtime.host.o -L../../bin/host -lhstreams_source


HOW TO RUN THE MATMUL  Block Matrix Multiplier REFERENCE CODE

1. Examine the runit.sh and ../common/setEnv.sh scripts.

You may need to edit it, depending on where you have the Intel
Composer XE compiler installed on your system.

Pay close attention to the setting for SINK_LD_LIBRARY_PATH, and
specifically the entries for /opt/mpss/ and the compiler.
There are multiple components:
  (a) mkl/lib/mic       : where to get the MKL libs for MIC side
  (b) compiler/lib/mic  : where to get the OpenMP libs for MIC side
  (c) /opt/mpss/3.4/sysroots/k1om-mpss-linux/usr/lib64 : where to get hstreams libs
  (d) /usr/lib64 : where to get host-side hstreams libs


2. Run the runit.sh script:

   $ . ./runit.sh 

It should nominally produce the following output:

~/my_hstreams/hstreams-release-3.4/bin/host ~/my_hstreams/hstreams-release-3.4/ref_code/matMult
~/my_hstreams/hstreams-release-3.4/bin/dev ~/my_hstreams/hstreams-release-3.4/bin/host ~/my_hstreams/hstreams-release-3.4/ref_code/matMult
~/my_hstreams/hstreams-release-3.4/bin/host ~/my_hstreams/hstreams-release-3.4/ref_code/matMult
~/my_hstreams/hstreams-release-3.4/ref_code/matMult
*** C(mxn) = A(mxk) x B(kxn)
blocksize = 500 x 500, m = 1000, k = 1500, n = 2000
mblocks = 2, kblocks = 3, nblocks = 4
Matrix in blocks A(2,3), B(3,4), C(2,4)
Perf: 32.934 Gflops/sec, Time= 182.185 msec
Perf: 220.094 Gflops/sec, Time= 27.261 msec
Perf: 253.969 Gflops/sec, Time= 23.625 msec
Perf: 256.094 Gflops/sec, Time= 23.429 msec
Size=, 1000, 1500, 1500, 2000, Pass=, 1, Iters=, 4, Max=, 256.09, GF/s, Avg_DGEMM=, 243.39, GFlop/s, StdDev=, 20.20, GFlop/s, 8.30, percent (Ignoring first iteration)
Computing result using host CPU...[If no MKLdgemm failures, then] Block Multiplication was successful.
MKL Host DGEMM Perf, 66.062, GFlops/sec, Time= 90.824 msec


HOW TO INTERPRET RESULTS OF THE MATMUL REFERENCE CODE

Note in the above, an estimated Gflops is computed.

To create a csv file to take to a spreadsheet program, do runs at different sizes saving output.txt file, then

    $ grep Size= output.txt > plot.csv 
