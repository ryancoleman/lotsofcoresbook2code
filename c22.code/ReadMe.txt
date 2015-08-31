========================================================================
    Gems_Code - build and execution instructions
========================================================================

Pre-requisites
--------------
1. Intel(R) Composer XE 2015 SP1 update 2 or later
2. GCC 4.4
3. Intel(R) Many Integrated Core (MIC) Platform Software Stack (MPSS)
   version 3.4 or later
4. Bash shell

To build and run on Xeon and Intel Xeon Phi type 'make'

Note: The build files assume the Xeon processor supports AVX instructions.
Readers with older CPUs will want to experiment with the '-xHost' Intel
compiler command line option (look for AVX in the makefiles). The -xHost
is a portable way to tell the compiler to use the fastest vector ISA
supported by the CPU.


Build and Makefile parameters
-----------------------------
Make utility file provided within a project is used to build the benchmark
for Intel(r) Xeon(r) and Intel(r) Xeon Phi(tm) hardware. In addition,
default invocation (w/o make parameters) invokes executables.

Makefile parameters:
- OMP=      select if vectorization shall be enabled through OpenMP* 4.0
            yes(default) - enable vectorization
			no           - generate scalar version of binaries
- PLAT=     select platform to build for.
            cpu(default) - build for Intel Xeon family processors
            knc          - build for Intel(r) Xeon Phi(tm) processor
- COMP=		select compiler family to be used, options:
            intel(default) - use Intel Composer XE
			gnu   - Use GNU C/C++ compiler


Benchmark execution
-------------------
Benchmark scripts:
- run_xeon.sh, this script executes performance evaluation on Intel Xeon processor
- run_mic.sh, this script executes performance evaluation on Intel Xeon Phi co-processor


/////////////////////////////////////////////////////////////////////////////
