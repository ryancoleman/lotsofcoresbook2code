This is the directory with the Original cpp_wilson_dslash 
code. 

Building the example code
========================

First Ensure Dependencies are built
then:

For Xeon Phi
------------

i) Build for Xeon Phi, by running ./build_mic.sh
ii) Run the timing on Xeon Phi by 
     1) Modifying ./run_test_mic.sh to suit your environment
     2) Launching the run with ./run_test_mic.sh

For Xeon 
---------

i) Build for Xeon, by running ./build.sh
ii) Run the timing on the Xeon by 
     1) Modifying ./run_test_xeon.sh to suit your environment
     2) Launching the run with ./run_test_xeon.sh

For Xeon with SSE Enabled
-------------------------
i) Build for Xeon with SSE enabled by running ./build_sse.sh
ii) Run the timing by 
     1) Modifying the ./run_test_xeon_sse.sh to suit your environmnet
     2) Launching the code with ./run_test_xeon_sse.sh


Notes:
======
*)  The build_xxx.sh scripts simply run configure; make; make check to generate the tests. The key differences are:
   a) build_mic.sh sources gets compiler flags from ../../env_mic.sh
      whereas the xeon targets use ../../env_avx.sh
   b) the Xeon SSE configuration has extra configure flags compared to the 
      regular Xeon configuration: --enable-sse2 --enable-sse3

*) The builds take place in the ./build_xxx subdirectories. 
   Executables are plased in the ./build_xxx/tests/

*) The only difference between run_test_xeon.sh and run_test_xeon_sse.sh is the location of the executable

*) run_test_mic.sh launches the job from the host using mpiexec.hydra. Your mileage may vary depending on your system

Contents:
=========

./build_mic.sh  -- Build Script for Xeon Phi
./build.sh -- Build Script for Xeon (same C code branch as for Xeon Phi)
./build_sse.sh --  Build Script for Xeon (SSE accelerated code branch)

./run_test_mic.sh       -- Run the timing on a Xeon Phi
./run_test_xeon.sh      -- Run the timing of the C-code on a Xeon 
./run_test_xeon_sse.sh  -- Run the timing of the SSE code branch on Xeon

./clean.sh -- Deletes everything in the build directories

cpp_wilson_dslash/ the source code
 +-- tests/   -- The testvol.h file for setting test volumes is in this 
 |                directory
 +-- lib/     -- the dispatch_scalar_openmp.cc and cpp_dslash_scalar_32bit.cc
 |                files live in this directory
 |
 +-- include/ -- cpp_dslash_scalar_32bit_c.h and cpp_dslash_scalar_32bit_sse.h 
                 files live in this direcotory

build_avx/    -- build.sh build the code in this direcotry for Xeon 
 +---tests/   -- the t_dslash and time_dslash code get built here

build_mic/    -- build_mic.sh builds the code in this directory for Xeon Phi
 +---tests/   -- the t_dslash and time_dslash code get built here

build_avx_sse/  -- build_sse.sh builds the code in this directory for Xeon
 +               with the SSE branch enabled
 +---tests/     -- the t_dslash and time_dslash code get built here

