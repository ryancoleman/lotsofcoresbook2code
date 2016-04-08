To build and run simply type 'make'

Note: Readers with older host CPUs will see a core dump at the very end
because some host benchmarks require AVX instructions. This only affects
the host openmp results.

Note: Readers who do not have their home directories automounted will
need to manually run the mic benchmarks (after typing 'make'):
   sh run_mic_nohome.sh bin.mic/benchmark_omp


Note:
=====

You can find the latest version on GitHub:
	https://github.com/noma/hexciton_benchmark

Please open an issue there or send an e-mail to ma.noack.pr@gmail.com if you
have any feedback on the benchmark.

Prepare:
========

Dependencies:
	- Intel C++ Compiler (the OpenCL examples can also be built with other
	  compilers, e.g. g++)
	- an OpenCL SDK
	- a Xeon Phi Coprocessor with MPSS installed (otherwise the MIC examples
	  will not build/run)
	+ everything listed in "thirdparty/Readme"
	
	For Plotting:
	- GNU R
	- pdfcrop
	
Follow the "Readme" in "thirdparty".

Build:
======

Use "make_omp.sh" and "make_ocl.sh" to build the OpenMP 4.0 and OpenCL
benchmarks. Have a look at the scripts, or run them with "-h" for further
details.

For the OpenCL benchmarks:
	./make_ocl.sh -c # Build the CPU version
	./make_ocl.sh -a # Build the Xeon Phi version
	./make_ocl.sh -g # Build the GPU version

For the OpenMP benchmarks:
	./make_omp.sh -c # Build the CPU version
	./make_omp.sh -a # Build the Xeon Phi version

Run:
====

OpenCL: 
	NOTE: Run from this directory, otherwise, the kernel sources won't be
	      found.
	Host (CL_DEVICE_TYPE_CPU): 
		bin/benchmark_ocl_host
	MIC (CL_DEVICE_TYPE_ACCELERATOR): 
		bin/benchmark_ocl_mic
	GPU (CL_DEVICE_TYPE_GPU): 
		bin/benchmark_ocl_gpu


OpenMP: 
	NOTE: If you do not have a shared home between host and MIC, you might have
	      to manually copy files to the MIC file system and manually start the
	      benchmark.
	Host:
		bin/benchmark_omp
	MIC (native):
		./run_mic.sh bin.mic/benchmark_omp

Evaluate:
=========

All benchmark results are written to standard output, all other messages to
standard error. Piping the output into a file will produce a tab-separated
easy-to-work-with table. Examples:

# Benchmark results are piped to results, meesages go to error.log:
bin/benchmark_ocl_mic > ocl_mic.data 2> error.log

# For better readability, pipe the output through column:
cat ocl_mic.data | column -t

# If you are only interested in a subset of the columns, use cut:
cut -f 1,2,4,5 ocl_mic.data | column -t

There are also GNU R scripts for generating PDF plots from the results:

# Plot benchmark results:
Rscript plot_results.r ocl_mic.data ocl_mic.pdf average

# Compare two benchmarks (from the same benchmark binary):
Rscript plot_compare.r ocl_mic_host1.data ocl_mic_host2.data ocl_host1_vs_host2.pdf average

APPENDIX:
=========

Assignment of plot labels to kernel source files.

Figure: runtime improvements per OpenCL optimisation
	"initial kernel"          "commutator_ocl_initial.cl"
	"first refactoring"	      "commutator_ocl_refactored.cl"
	"AoSoA naive"             "commutator_ocl_aosoa_naive.cl"
	"AoSoA 2D-NDRange"        "commutator_ocl_aosoa.cl"
	"compile-time constants"  "commutator_ocl_aosoa_constantscl"
	"manual vectorization"	  "commutator_ocl_manual_aosoa_constants_direct_perm.cl" (CPU)
                          "commutator_ocl_manual_aosoa_constants_direct_prefetch.cl" (MIC)
	
Figure: OpenCL vectorization approaches
	"AoSoA 2D-NDRange"                 "commutator_ocl_aosoa.cl"
		                               "commutator_ocl_manual_aosoa.cl"
	"compile-time constants"           "commutator_ocl_aosoa_constants.cl"
		                               "commutator_ocl_manual_aosoa_constants.cl"
	"direct accumulation"              "commutator_ocl_aosoa_direct.cl"
		                               "commutator_ocl_manual_aosoa_direct.cl"
	"constants + direct write"         "commutator_ocl_aosoa_constants_direct.cl"
		                               "commutator_ocl_manual_aosoa_constants_direct.cl"
	"constants, direct, permuted k,j"  "commutator_ocl_aosoa_constants_direct_perm.cl"
		                               "commutator_ocl_manual_aosoa_constants_direct_perm.cl"

Figure: performance portability across different devices
	"fastest kernel on Sandy Bridge"  "commutator_ocl_aosoa"
	"fastest kernel on Haswell"       "commutator_ocl_manual_aosoa_constants_direct_perm"
	"fastest kernel on MIC"           "commutator_ocl_manual_aosoa_constants_direct_prefetch"
	"fastest kernel on GPGPU"         "commutator_ocl_gpu_final"

Figure: OpenMP kernel runtimes
	"AoSoA base kernel, auto"                                 "commutator_omp_aosoa.cpp"
	"AoSoA base kernel, manual"                               "commutator_omp_manual_aosoa.cpp"
	"constants + direct, auto"                                "commutator_omp_aosoa_constants_direct.cpp"
	"constants + direct, manual"                              "commutator_omp_manual_aosoa_constants_direct.cpp"
	"constants + direct\n+ loop 2->3, auto"                   "commutator_omp_aosoa_constants_direct_perm2to3.cpp"
	"constants + direct + loop 4->5\n+ unroll hints, manual"  "commutator_omp_manual_aosoa_constants_direct_perm_unrollhints.cpp"

