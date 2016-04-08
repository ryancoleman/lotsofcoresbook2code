SHOC Alpha v1.1.4a-mic
======================

Alpha Status
------------

This release is an alpha.  It is expected to have bugs and has been tested on
a limited number of platforms compared to the main distribution.  The
supported platform is CentOS 6.2 \ Intel Compiler v13.0.1

Please contact shoc-help@elist.ornl.gov with any questions concerning this
release.

Please see the release notes below for other limitations and known issues.

Building the Benchmarks
-----------------------

Type ```make``` inside this directory. This assumes the intel toolchain 
(a recent version of ```icc```) is in your current path. Your paths
should also include MKL, which is used in some of the benchmarks.

Running the Benchmarks
----------------------

How to run all the benchmarks:

1) Use the provided run script in ```/bin```

How to run an individual benchmark:

2) Set the appropriate environment variables (see run scripts)

For example, to run MD workload
```
    $ export MIC_ENV_PREFIX=MIC
    $ export MIC_USE_2MB_BUFFERS=32K
    $ export MIC_KMP_AFFINITY=granularity=fine,balanced
    $ export MIC_BUFFERSIZE=128M
    $ export nc=61
    $ export MIC_OMP_NUM_THREADS=4*$nc-4
    $ ./MD -s 4
```
For another example, to run FFT workload
```
    $ export MIC_ENV_PREFIX=MIC 
    $ export MIC_USE_2MB_BUFFERS=32K 
    $ export MIC_KMP_AFFINITY=granularity=fine,balanced 
    $ export MIC_BUFFERSIZE=128M 
    $ export MIC_MKL_DYNAMIC=false
    $ nc=61
    $ export MIC_OMP_NUM_THREADS=4*$nc-4
    $ ./FFT -s 4
```

The value for ```-s``` is a number between 1 and four that corresponds to 
problem size from smallest to largest (consistent with the CUDA and 
OpenCL versions).

Release Notes
-------------

  * This is a preview release and does not support MPI or more than one MIC card
per node.
  * The current benchmarks have not yet been integrated into the main SHOC
infrastructure.  As such, a simplified runscript has been provided in the /bin
directory.
  * MIC versions of the benchmarks have can have (much) longer build\run times
than CUDA/OCL. In particular, the compile time for S3D and the runtime for
DeviceMemory are quite long. 

Known Issues
------------

  * DeviceMemory is inconsistent with CUDA/OpenCL.  Its tests are structured to
be specific for MIC hardware.
  * Current release uses gettimeofday() exclusively. In most cases this shouldn't
matter.
  * Sort does not currently complete all radix passes.

