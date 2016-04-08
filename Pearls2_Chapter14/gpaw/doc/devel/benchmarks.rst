.. _benchmarks:

==========
Benchmarks
==========

.. _memory_bandwidth:

Memory benchmark
================

Goal
----

It is known that GPAW puts a heavy load on the RAM memory subsystem.
This benchmark will test system's
`memory bandwidth <http://en.wikipedia.org/wiki/Memory_bandwidth>`_.

Prerequisites
-------------

This benchmark requires approximately 1.5 GB of RAM memory per core.
The amount of disk space required is minimal.

The following packages are required (names given for RHEL 5 system):

 - python, python-devel
 - numpy
 - python-matplotlib
 - openmpi, openmpi-devel
 - bash
 - `campos-gpaw <https://wiki.fysik.dtu.dk/gpaw/install/installationguide.html>`_
 - `campos-ase3 <https://wiki.fysik.dtu.dk/ase/download.html>`_

Please refer to :ref:`platforms_and_architectures` for hints on
installing GPAW on different platforms.

Results
-------

Multiple instances of the GPAW code are executed in serial
using OpenMPI in order to benchmark a number of processes that ranges from
1, through integer powers of 2 and up to the total number of CPU cores
(NCORES - number of cores available on the test machine).

The benchmark result is the average execution time in seconds when running
1, 2, up to NCORES processes, respectively, on the test machine.
The scaling of the execution time with the number of processes is part of
the benchmark result.

On `Non-Uniform Memory Access <http://en.wikipedia.org/wiki/Non-Uniform_Memory_Access>`_ machines,
in order to get reproducible results,
it is important to find the (fastest) combination of processor/memory.
See section :ref:`opteron_285` for an example of performance degradation
depending on the numa policy.

To see which CPU corresponds to which node examine
``/sys/devices/system/node/node*``.  It is assumed the following
mapping:

- dual-socket dual-core machine::

   numactl --membind=0 --physcpubind=0
   numactl --membind=0 --physcpubind=1
   numactl --membind=1 --physcpubind=2
   numactl --membind=1 --physcpubind=3

- dual-socket quad-core machine::

   numactl --membind=0 --physcpubind=0
   numactl --membind=0 --physcpubind=1
   numactl --membind=0 --physcpubind=2
   numactl --membind=0 --physcpubind=3
   numactl --membind=1 --physcpubind=4
   numactl --membind=1 --physcpubind=5
   numactl --membind=1 --physcpubind=6
   numactl --membind=1 --physcpubind=7

**Note** that the mapping above assigns first all the cores belonging to the
first memory node, what is **non**-optimal for runs with the number
of processes larger than 1 and smaller than NCORES.
See :ref:`opteron_285` section for example of this behaviour.

For description of numa see `NUMACTL <https://computing.llnl.gov/LCdocs/chaos/index.jsp?show=s5.2.2>`_
and `NUMA support <http://lwn.net/Articles/254445/>`_.

Getting the results
-------------------

Please perform the following steps:

 - make sure that no other resources consuming processes are running,
 - set (as root) ulimit's cpu time to 5 hours::

    ulimit -t 18000

 - use the following commands to setup the benchmark::

    bash
    export NCORES=8 # default; set this variable to the number of cores in your machine
    export CORES_PER_SOCKET=4 # default; set this variable to the number of cores per socket
    export MACHINE=TEST # default; optional: set this to the name of you machine
    mkdir /tmp/benchmark.$$; cd /tmp/benchmark.*
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/memory_bandwidth/H2Al110.py
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/memory_bandwidth/prepare.sh
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/memory_bandwidth/taskit.BINDING.one.node
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/memory_bandwidth/run.sh
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/memory_bandwidth/run_numactl.sh
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/memory_bandwidth/memory_bandwidth.py
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/memory_bandwidth/twiny.py
    chmod u+x taskit.BINDING.one.node
    sh prepare.sh

 - run with (it takes 6-10 hours with NCORES=8)::

    cd $MACHINE; nohup sh ../run.sh 2>&1 | tee $MACHINE.log&

   **Warning**: on numa-enabled machines use::

    cd $MACHINE; nohup sh ../run_numactl.sh 2>&1 | tee $MACHINE.log&

 - analyse the results::

    python ../memory_bandwidth.py

 - to estimate performance run the benchmark on the maximal number the cores only::

    export NCORES=8
    export CORES_PER_SOCKET=4
    export MACHINE=TEST
    export STARTCORES=${NCORES}
    cd $MACHINE; nohup sh ../run_numactl.sh 2>&1 | tee $MACHINE.log&

Benchmarked systems
-------------------

.. _best_performance:

Best performance
++++++++++++++++

The best performance estimate has been obtained on the following systems
with the following configuration of GPAW tested on **production runs**:
compiler:blas/lapack/(numpy:dotblas/lapack):

- GPAW **0.7.6383** (28 SCF steps):

  - Xeon_X5570_: 329.0 s (11.8 s/step) - gcc43:goto2-1.13/acml-4.4.0/(numpy:default/acml-4.0.1),
    date: May 08 2010, kernel 2.6.18-128.7.1.el5, BIOS HP O33 02/04/2010.

  - opteron_285_: 659.9 s (23.6 s/step) - gcc43:goto-1.26/acml-4.4.0/(numpy:default/acml-4.0.1)*,
    date: May 08 2010, kernel 2.6.18-164.el5, BIOS IBM 1.35 07/18/2007.

- GPAW **0.6.5147** (30 SCF steps):

  - Xeon_X5570_: 345.1 s (11.5 s/step) - gcc43:goto2-1.13/acml-4.4.0/(numpy:default/acml-4.0.1),
    date: May 08 2010, kernel 2.6.18-128.7.1.el5, BIOS HP O33 02/04/2010.

  - Xeon_X5667_: 509.8 s (17.0 s/step) - gcc43:acml-4.3.0/acml-4.3.0/(numpy:default/acml-4.0.1),
    date: May 08 2010, kernel 2.6.18-164.15.1.el5, BIOS HP 0.34 03/31/2010.

  - opteron_285_: 674.9 s (22.5 s/step) - gcc43:goto-1.26/acml-4.3.0/(numpy:default/acml-4.0.1)*,
    date: May 08 2010, kernel 2.6.18-164.el5, BIOS IBM 1.35 07/18/2007.

See the above links for the detailed results.

.. _opteron_285:

Dual-socket dual Core AMD Opteron(tm) Processor 285/2.6 GHz/2 GB RAM per core EL5
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- memory bandwidth:

  - date: May 08 2010, kernel 2.6.18-164.el5, BIOS IBM 1.35 07/18/2007.
    Performed with gcc43/goto-1.26/acml-4.4.0, GPAW **0.7.6383** (28 SCF steps),
    numpy *1.3.0* compiled with gcc/default(no dotblas)/acml-4.0.1:

    - run with assumed numactl mapping for a dual-socket dual-core machine::

       export NCORES=4
       export CORES_PER_SOCKET=4
       export MACHINE=gcc43.numactl
       export STARTCORES=${NCORES}
       cd $MACHINE; nohup sh ../run_numactl.sh 2>&1 | tee $MACHINE.log&

     results in::

       No. of processes 1: time [sec]: avg 564.4, stddev 0.6, min 563.4, max 565.1
       No. of processes 2: time [sec]: avg 658.0, stddev 3.6, min 653.0, max 662.9
       No. of processes 4: time [sec]: avg 659.9, stddev 3.8, min 654.4, max 666.1

  - date: May 08 2010, kernel 2.6.18-164.el5, BIOS IBM 1.35 07/18/2007.
    Performed with gcc43/goto-1.26/acml-4.3.0, GPAW **0.6.5147** (30 SCF steps),
    numpy *1.3.0* compiled with gcc/default(no dotblas)/acml-4.0.1:

    - run with assumed numactl mapping for a dual-socket dual-core machine::

       export NCORES=4
       export CORES_PER_SOCKET=4
       export MACHINE=gcc43.numactl
       export STARTCORES=${NCORES}
       cd $MACHINE; nohup sh ../run_numactl.sh 2>&1 | tee $MACHINE.log&

     results in::

       No. of processes 1: time [sec]: avg 586.6, stddev 0.7, min 585.7, max 587.7
       No. of processes 2: time [sec]: avg 673.9, stddev 3.4, min 669.4, max 678.5
       No. of processes 4: time [sec]: avg 674.9, stddev 3.2, min 671.1, max 681.9

  - date: ??, kernel ??, BIOS IBM ??.
    Performed with gcc43/goto-1.26/acml-4.2.0, GPAW **0.6.3862** (35 SCF steps),
    numpy *1.3.0* compiled with gcc/blas-3.0-37/lapack-3.0-37:

    - run with default numa::

       export NCORES=4
       export CORES_PER_SOCKET=2
       export MACHINE=gcc43
       export STARTCORES=${NCORES}
       cd $MACHINE; nohup sh ../run.sh 2>&1 | tee $MACHINE.log&

     results in::

      No. of processes 1: time [sec]: avg 716.1, stddev 3.7, min 710.8, max 719.6
      No. of processes 2: time [sec]: avg 726.9, stddev 7.2, min 718.2, max 735.0
      No. of processes 4: time [sec]: avg 898.6, stddev 7.5, min 890.5, max 914.1

    - run with assumed numactl mapping for a dual-socket dual-core machine::

       export NCORES=4
       export CORES_PER_SOCKET=2
       export MACHINE=gcc43.numactl
       export STARTCORES=${NCORES}
       cd $MACHINE; nohup sh ../run_numactl.sh 2>&1 | tee $MACHINE.log&

     results in::

      No. of processes 1: time [sec]: avg 717.5, stddev 0.8, min 716.0, max 718.1
      No. of processes 2: time [sec]: avg 884.7, stddev 7.7, min 873.4, max 897.1
      No. of processes 4: time [sec]: avg 894.3, stddev 15.4, min 874.9, max 913.9

    **Note** the performance degradation in the case of numactl for two cores,
    compared to a "default" run. The degradation of ~25% between 1 core and the maximal number
    of cores (4) is typical for this generation of AMD systems.

- performance estimate (average time of the memory_bandwidth_ benchmark on the maximal number of cores):

  - GPAW **0.6.5147** (30 SCF steps) was used.
    Standard deviations are found below 15 sec. "**N/A**" denotes the fact that libraries are not available,
    "**-**" that tests were not performed.

    =================================================== ========= ========= ========= ============
    blas/lapack/(numpy:dotblas/lapack): compiler        gcc       gcc43     icc 11.0  open64 4.2.3
    =================================================== ========= ========= ========= ============
    acml-4.4.0/acml-4.4.0/(default/acml-4.0.1)*         N/A        716.4    --         689.3
    acml-4.4.0/acml-4.4.0/(blas-3.0-37/lapack-3.0-37)   N/A        --       --         669.0
    acml-4.3.0/acml-4.3.0/(default/acml-4.0.1)*         N/A        713.5    --        --
    acml-4.3.0/acml-4.3.0/(blas-3.0-37/lapack-3.0-37)   N/A        699.7    --        --
    acml-4.0.1/acml-4.0.1/(default/acml-4.0.1)*          715.4    N/A       --        --
    blas-3.0-37/lapack-3.0-37/(default/acml-4.0.1)*     1151.6 F  1146.3 F  --        --
    goto2-1.13/acml-4.4.0/(default/acml-4.0.1)*         N/A        680.4        --     652.9
    goto2-1.13/acml-4.4.0/(blas-3.0-37/lapack-3.0-37)   N/A        699.6        --     669.0
    goto-1.26/acml-4.4.0/(default/acml-4.0.1)*          N/A        680.4        --     651.1
    goto-1.26/acml-4.3.0/(default/acml-4.0.1)*          N/A        674.9    --        --
    goto-1.26/acml-4.3.0/(blas-3.0-37/lapack-3.0-37)    N/A        693.2    --        --
    atlas-3.8.3/atlas-3.8.3/(default/acml-4.0.1)*       --        FAIL      --        --
    =================================================== ========= ========= ========= ============

    **Note**: the numpy version marked by \* (star) denotes that the ``_dotblas.so``
    module was **NOT** built and the given lapack used.

    **Warning**: fields marked by **F** denote a failure in the GPAW's test suite.
    Fields marked by **FAIL** denote a failure in the memory_bandwidth_ benchmark.
    Errors were reported when using different blas/lapack in GPAW and NUMPY!

    ============================== =============================================
    compiler                       options                     
    ============================== =============================================
    gcc 4.1.2 20080704             -O3 -funroll-all-loops -std=c99
    gcc43 4.3.2 20081007           -O3 -funroll-all-loops -std=c99
    icc 11.0 083                   -xHOST -O3 -ipo -no-prec-div -static -std=c99
    open64 4.2.3                   -O3 -std=c99 -fPIC
    ============================== =============================================

  - GPAW **0.6.3862** (35 SCF steps) was used, numpy *1.3.0* compiled with gcc/goto-1.26/acml-4.0.1.
    Standard deviations are found below 15 sec. "**N/A**" denotes the fact that libraries are not available,
    "**-**" that tests were not performed.

    ============================= ======= ======= ======= ======= ======= =======
    blas/lapack : compiler        gcc     gcc43   amd4.2  pathcc  icc     pgcc
    ============================= ======= ======= ======= ======= ======= =======
    acml-4.2.0/acml-4.2.0         N/A      991.74  985.83  980.75 1020.66 1082.64
    acml-4.1.0/acml-4.1.0         N/A     --      --       978.58 --      --     
    acml-4.0.1/acml-4.0.1          991.95 N/A     N/A      984.23 --      --     
    blas-3.0-37/lapack-3.0-37     1494.63 1495.52 --      --      --      --     
    goto-1.26/acml-4.2.0          N/A      889.22  886.43 879.28  FAIL    FAIL
    goto-1.26/acml-4.2.0 PGO      --       886.47 --      --      --      --     
    goto-1.26/acml-4.0.1           888.79 N/A     N/A     --      --      --     
    atlas-3.8.3/acml-4.2.0        --       931.41 --      --      --      --     
    atlas-3.8.3/lapack-3.2.1      --       927.71 --      --      --      --     
    mkl-10.1.2.024/mkl-10.1.2.024 --      1012.64 --      1030.06 --      --     
    ============================= ======= ======= ======= ======= ======= =======

    **Note**: the PGO entry refers to :ref:`PGO` driven using the benchmark.

    **Warning**: fields marked by **FAIL** denote a failure in the memory_bandwidth_ benchmark.
    Errors were reported when using different blas/lapack in GPAW and NUMPY!

    ============================== =============================================
    compiler                       options                     
    ============================== =============================================
    gcc 4.1.2 20080704             -O3 -funroll-all-loops -std=c99
    gcc43 4.3.2 20081007           -O3 -funroll-all-loops -std=c99
    gcc 4.2.0-amd-barcelona-rhel4  -O3 -funroll-all-loops -std=c99
    pathcc Version 3.2 2008-06-16  -O3 -OPT:Ofast -ffast-math -std=c99
    icc 11.0 083                   -xHOST -O3 -ipo -no-prec-div -static -std=c99
    pgcc 8.0-6                     -fast -Mipa=fast,inline -Msmartalloc
    ============================== =============================================

    **Note**: that using wrong numa policy (in some situations also the **default** numa policy)
    results in serious performance degradation, and non-reproducible results.
    Example below is given for gcc43/goto-1.26/acml-4.2.0 (**A**),
    and gcc43/mkl-10.1.2.024/mkl-10.1.2.024 (**B**).

    ===================== ====================================== ======================================
    MP pairs (see below)  A Runtime [sec]                        B Runtime [sec]                       
    ===================== ====================================== ======================================
    00 01 12 13           avg 889.22, stddev 7.61, max 902.53    avg 1012.64, stddev 11.65, max 1032.98
    default (not set)     avg 892.22, stddev 12.54, max 915.96   avg 1047.2, stddev 51.8, max 1171.5
    00 11 02 13           avg 953.39, stddev 81.57, max 1069.16  avg 1081.78, stddev 92.67, max 1204.43
    10 11 02 03           avg 1330.88, stddev 11.75, max 1351.37 avg 1504.35, stddev 8.89, max 1527.54
    00 01 02 03           avg 1549.29, stddev 59.61, max 1645.92 avg 1736.57, stddev 77.87, max 1849.49
    ===================== ====================================== ======================================

    **Note**: "MP pairs" denote pairs of M and P used for `numactl --membind=M --physcpubind=P`
    for ranks 0, 1, 2, 3, respectively.
    In this case **A** the **default** numa policy does not result in performance degradation.

.. _Xeon_X5570:

Dual-socket quad Core 64-bit Intel Nehalem Xeon X5570 quad-core 2.93 GHz/3 GB RAM per core EL5
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- memory bandwidth:

  - date: May 08 2010, kernel 2.6.18-128.7.1.el5, BIOS HP O33 02/04/2010.
    Performed with gcc43/goto2-1.13/acml-4.4.0, GPAW **0.7.6383** (28 SCF steps),
    numpy *1.3.0* compiled with gcc/default(no dotblas)/acml-4.0.1:

    - run with assumed numactl mapping for a dual-socket quad-core machine::

       export NCORES=8
       export CORES_PER_SOCKET=4
       export MACHINE=gcc43.numactl
       export STARTCORES=${NCORES}
       cd $MACHINE; nohup sh ../run_numactl.sh 2>&1 | tee $MACHINE.log&

     results in::

      No. of processes 1: time [sec]: avg 297.4, stddev 0.3, min 296.8, max 297.7
      No. of processes 2: time [sec]: avg 307.0, stddev 0.9, min 305.8, max 308.6
      No. of processes 4: time [sec]: avg 327.9, stddev 0.9, min 326.5, max 329.6
      No. of processes 6: time [sec]: avg 321.7, stddev 10.3, min 306.3, max 330.7
      No. of processes 8: time [sec]: avg 329.0, stddev 1.5, min 326.9, max 332.5

  - date: May 08 2010, kernel 2.6.18-128.7.1.el5, BIOS HP O33 02/04/2010.
    Performed with gcc43/goto2-1.13/acml-4.4.0, GPAW **0.6.5147** (30 SCF steps),
    numpy *1.3.0* compiled with gcc/default(no dotblas)/acml-4.0.1:

    - run with assumed numactl mapping for a dual-socket quad-core machine::

       export NCORES=8
       export CORES_PER_SOCKET=4
       export MACHINE=gcc43.numactl
       export STARTCORES=${NCORES}
       cd $MACHINE; nohup sh ../run_numactl.sh 2>&1 | tee $MACHINE.log&

     results in::

      No. of processes 1: time [sec]: avg 313.2, stddev 0.2, min 313.0, max 313.6
      No. of processes 2: time [sec]: avg 322.9, stddev 1.2, min 321.5, max 324.9
      No. of processes 4: time [sec]: avg 344.1, stddev 0.8, min 342.5, max 345.7
      No. of processes 6: time [sec]: avg 337.5, stddev 10.1, min 322.8, max 347.8
      No. of processes 8: time [sec]: avg 345.1, stddev 1.5, min 343.1, max 348.9

- performance estimate (average time of the memory_bandwidth_ benchmark on the maximal number of cores):

  - GPAW **0.6.5147** (30 SCF steps) was used.
    Standard deviations are found below 15 sec. "**N/A**" denotes the fact that libraries are not available,
    "**-**" that tests were not performed.

    ============================================================= ========= ========= ========= ========= ============
    blas/lapack/(numpy:dotblas/lapack): compiler                  gcc       gcc43     icc 11.0  icc 11.1  open64 4.2.3 
    ============================================================= ========= ========= ========= ========= ============
    acml-4.4.0/acml-4.4.0/(default/acml-4.0.1)*                   N/A        436.6     399.2 F   400.0 F   418.5
    acml-4.4.0/acml-4.4.0/(blas-3.0-37/lapack-3.0-37)             N/A        355.5     326.7 F   326.0 F   347.4
    acml-4.3.0/acml-4.3.0/(default/acml-4.0.1)*                   N/A        435.9    --        --        --
    acml-4.3.0/acml-4.3.0/(blas-3.0-37/lapack-3.0-37)             N/A        364.8    --        --        --
    acml-4.3.0/acml-4.3.0/(mkl-10.1.3.027/mkl-10.1.3.027)         N/A        363.4    --        --        --
    acml-4.0.1/acml-4.0.1/(default/acml-4.0.1)*                    443.5    N/A       --        --        --
    blas-3.0-37/lapack-3.0-37/(default/acml-4.0.1)*                529.7  F  531.2 F  --        --        --
    goto2-1.13/acml-4.4.0/(default/acml-4.0.1)*                   N/A        345.1        --        --     326.6
    goto2-1.13/acml-4.4.0/(blas-3.0-37/lapack-3.0-37)             N/A        351.1        --        --     333.3
    goto-1.26/acml-4.3.0/(default/acml-4.0.1)*                    N/A       N/A       N/A       N/A       N/A
    atlas-3.8.3/atlas-3.8.3/(default/acml-4.0.1)*                 --         380.0 F  --        --        --
    mkl-10.1.3.027/mkl-10.1.3.027/(default/acml-4.0.1)*           --         352.3     318.4 F  --        --
    mkl-10.1.3.027/mkl-10.1.3.027/(mkl-10.1.3.027/mkl-10.1.3.027) --         382.9     332.4 F  --        --
    mkl-10.1.3.027/mkl-10.1.3.027/(blas-3.0-37/lapack-3.0-37)     --         358.0     326.5 F  --        --
    ============================================================= ========= ========= ========= ========= ============

    **Note**: the numpy version marked by \* (star) denotes that the ``_dotblas.so``
    module was **NOT** built and the given lapack used.

    **Warning**: fields marked by **F** denote a failure in the GPAW's test suite.
    Fields marked by **FAIL** denote a failure in the memory_bandwidth_ benchmark.
    Errors were reported when using different blas/lapack in GPAW and NUMPY!

    ============================== =============================================
    compiler                       options                     
    ============================== =============================================
    gcc 4.1.2 20080704             -O3 -funroll-all-loops -std=c99
    gcc43 4.3.2 20081007           -O3 -funroll-all-loops -std=c99
    icc 11.0 083                   -xHOST -O3 -ipo -no-prec-div -static -std=c99
    icc 11.1 072                   -xHOST -O3 -ipo -no-prec-div -static -std=c99
    open64 4.2.3                   -O3 -std=c99 -fPIC
    ============================== =============================================

.. _Xeon_X5667:

Dual-socket quad Core 64-bit Intel Westmere Xeon X5667 quad-core 3.07 GHz/3 GB RAM per core EL5
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

**Note**: the benchmark was performred with a pre-release system of CPU and beta-version of BIOS.
The performance numbers may not reflect the future production systems.

- memory bandwidth:

  - date: May 08 2010, kernel 2.6.18-164.15.1.el5, BIOS HP 0.34 03/31/2010.
    Performed with gcc43/acml-4.3.0/acml-4.3.0, GPAW **0.6.5147** (30 SCF steps),
    numpy *1.3.0* compiled with gcc/default(no dotblas)/acml-4.0.1:

    - run with assumed numactl mapping for a dual-socket quad-core machine::

       export NCORES=8
       export CORES_PER_SOCKET=4
       export MACHINE=gcc43.numactl
       export STARTCORES=${NCORES}
       cd $MACHINE; nohup sh ../run_numactl.sh 2>&1 | tee $MACHINE.log&

     results in::

       No. of processes 1: time [sec]: avg 423.3, stddev 0.9, min 422.2, max 424.8
       No. of processes 2: time [sec]: avg 452.7, stddev 0.5, min 451.9, max 453.5
       No. of processes 4: time [sec]: avg 481.0, stddev 1.7, min 479.0, max 484.5
       No. of processes 6: time [sec]: avg 483.1, stddev 13.6, min 462.3, max 497.3
       No. of processes 8: time [sec]: avg 509.8, stddev 2.5, min 506.7, max 517.1

- performance estimate (average time of the memory_bandwidth_ benchmark on the maximal number of cores):

  - GPAW **0.6.5147** (30 SCF steps) was used.
    Standard deviations are found below 15 sec. "**N/A**" denotes the fact that libraries are not available,
    "**-**" that tests were not performed.

    ============================================================= ========= ========= ========= ========= =========
    blas/lapack/(numpy:dotblas/lapack): compiler                  gcc       gcc43     gcc44     icc 11.0  icc 11.1
    ============================================================= ========= ========= ========= ========= =========
    acml-4.4.0/acml-4.4.0/(default/acml-4.0.1)*                   N/A       --         511.0     469.4 F   469.8 F
    acml-4.4.0/acml-4.4.0/(blas-3.0-37/lapack-3.0-37)             N/A        454.2     --        418.6 F   419.0 F
    acml-4.3.0/acml-4.3.0/(default/acml-4.0.1)*                   N/A        509.8     510.8    --        --
    acml-4.3.0/acml-4.3.0/(blas-3.0-37/lapack-3.0-37)             N/A        454.3     453.3    --        --
    acml-4.3.0/acml-4.3.0/(mkl-10.1.3.027/mkl-10.1.3.027)         N/A        452.9    --        --        --
    acml-4.0.1/acml-4.0.1/(default/acml-4.0.1)*                    508.6    N/A       N/A       --        --
    blas-3.0-37/lapack-3.0-37/(default/acml-4.0.1)*               --        --        --        --        --
    goto-1.26/acml-4.3.0/(default/acml-4.0.1)*                    N/A       N/A       N/A       N/A       N/A
    atlas-3.8.3/atlas-3.8.3/(default/acml-4.0.1)*                 --        --        --        --        --
    mkl-10.1.3.027/mkl-10.1.3.027/(default/acml-4.0.1)*           --         429.3    --        --        --
    mkl-10.1.3.027/mkl-10.1.3.027/(mkl-10.1.3.027/mkl-10.1.3.027) --         440.9    --         405.6 F  --
    mkl-10.2.1.017/mkl-10.2.1.017/(mkl-10.1.3.027/mkl-10.1.3.027) --         439.6 F  --        --        --
    mkl-10.2.4.032/mkl-10.2.4.032/(mkl-10.1.3.027/mkl-10.1.3.027) --         442.5 F   438.3 F  --        --
    mkl-10.1.3.027/mkl-10.1.3.027/(blas-3.0-37/lapack-3.0-37)     --         440.4    --        --        --
    ============================================================= ========= ========= ========= ========= =========

    **Note**: the numpy version marked by \* (star) denotes that the ``_dotblas.so``
    module was **NOT** built and the given lapack used.

    **Warning**: fields marked by **F** denote a failure in the GPAW's test suite.
    Fields marked by **FAIL** denote a failure in the memory_bandwidth_ benchmark.
    Errors were reported when using different blas/lapack in GPAW and NUMPY!

    ============================== =============================================
    compiler                       options                     
    ============================== =============================================
    gcc 4.1.2 20080704             -O3 -funroll-all-loops -std=c99
    gcc43 4.3.2 20081007           -O3 -funroll-all-loops -std=c99
    gcc44 4.4.0 20090514           -O3 -funroll-all-loops -std=c99
    icc 11.0 083                   -xHOST -O3 -ipo -no-prec-div -static -std=c99
    icc 11.1 072                   -xHOST -O3 -ipo -no-prec-div -static -std=c99
    ============================== =============================================

Strong scaling benchmarks
=========================

Goal
----

Fix the problem size, vary the number of processors, and measure the speedup.

1) Medium size system
+++++++++++++++++++++

The system used in this benchmark is of medium size, as for the year 2008,
and consists of 256 water molecules in a box of ~20**3 Angstrom**3,
2048 electrons, and 1056 bands, and 112**3 grid points (grid spacing of ~0.18).
LCAO initialization stage is performed, then 3 SCF steps with a constant
potential and 2 full SCF steps.
All the stages are timed separately, due to their different scaling.

**Note** that the size of the system can be changed easily by modifying
just one variable in :svn:`~doc/devel/256H2O/b256H2O.py`::

  r = [2, 2, 2]

Prerequisites
+++++++++++++

This benchmark requires approximately 2 GB of RAM memory per core and at least 16 cores.
The amount of disk space required is minimal.

The following packages are required (names given for Fedora Core 10 system):

 - python, python-devel
 - numpy
 - python-matplotlib
 - openmpi, openmpi-devel
 - blacs, scalapack
 - bash
 - `campos-gpaw <https://wiki.fysik.dtu.dk/gpaw/install/installationguide.html>`_
 - `campos-ase3 <https://wiki.fysik.dtu.dk/ase/download.html>`_

**Note** that GPAW has to built with ScaLAPACK enabled -
please refer to :ref:`platforms_and_architectures` for hints on
installing GPAW on different platforms.

Results
+++++++

GPAW code is executed in parallel in order to benchmark a number of processes that ranges from 16,
through integer powers of 2 up to 128.

The number of bands (1056) and cores are chosen to make comparisons
of different band parallelizations (:ref:`band_parallelization`) possible.

**Note**: to achive optimal performance diagonalization steps are performed
on `4x4` blacs grid with block size of `64` specified by adding ``--gpaw=blacs=1 --sl_diagonalize=4,4,64`` options.

**Note** also that a default domain decomposition is appplied, and different
results can be obtained by tuning ``--domain-decomposition`` argument
to your platform (see :ref:`submit_tool_on_niflheim`).

**Note**: the ``--gpaw=usenewlfc=1`` option is required to skip the calculation of forces
and decrease **memory** usage.

The results of the benchmark is scaling of execution time of different stages
of GPAW run with the number of processes (CPU cores).

Getting the results
+++++++++++++++++++

Please perform the following steps:

 - use the following commands to setup the benchmark::

    bash
    mkdir 256H2O; cd 256H2O
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/256H2O/b256H2O.py
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/256H2O/akka.sh
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/256H2O/surveyor.sh
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/256H2O/prepare.sh
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/256H2O/scaling.py
    # set the prefix directory: results will be in $PATTERN_*_
    export PATTERN=b256H2O_112_04x04m64.grid
    sh prepare.sh

   **Warning**: the choice of the directory names is not free in the sense that
   the number of processes has to come at the end of directory name,
   and be delimited by two underscores.

 - run with, for example:

    - on akka::

       cd ${PATTERN}_00016_; qsub -l nodes=2:8 ../akka.sh; cd ..
       cd ${PATTERN}_00032_; qsub -l nodes=4:8 ../akka.sh; cd ..
       cd ${PATTERN}_00064_; qsub -l nodes=8:8 ../akka.sh; cd ..
       cd ${PATTERN}_00128_; qsub -l nodes=16:8 ../akka.sh; cd ..

   **Warning**: on Linux clusters it s desirable to repeat these runs 2-3 times
   to make sure that they give reproducible time. Even with this procedure obtained
   runtimes may show up to 5% precision.

 - analyse the results::

    python scaling.py -v --dir=. --pattern="b256H2O_112_04x04m64.grid_*_" b256H2O

   Niflheim results:

   - opteron (IBM eServer x3455: Opteron 2218 dual-core 2.60 GHz CPUs) nodes (infiniband):
     performed on EL4 with gcc/acml-4.0.1/acml-4.0.1, GPAW **0.6.5092**,
     numpy *1.0.3* compiled with gcc/blas-3.0-25/lapack-3.0-25 (with dotblas); no ScaLAPACK used::

       # p - processes, p0 - reference processes, t - time [sec], s - speedup, e - efficiency
       # GPAW version 6.5092: stages: 1 - initialization, 2 - fixdensity, 3 - SCF, 4 - forces, 5 - total
       # p     p/p0   t1      s1      e1    t2      s2      e2    t3      s3      e3    t4      s4      e4    t5      s5      e5
            16   1.00   201.5    16.0  1.00   778.5    16.0  1.00   533.0    16.0  1.00     0.0     0.0  0.00  1513.0    16.0  1.00
            32   2.00   113.5    28.4  0.89   391.5    31.8  0.99   267.0    31.9  1.00     0.0     0.0  0.00   772.0    31.4  0.98
            64   4.00    69.0    46.7  0.73   204.0    61.1  0.95   139.0    61.4  0.96     0.0     0.0  0.00   412.0    58.8  0.92

   - opteron (IBM eServer x3455: Opteron 2218 dual-core 2.60 GHz CPUs) nodes (ethernet):
     performed on EL5 with gcc43/goto-1.26/acml-4.3.0, GPAW **0.6.5092**,
     numpy *1.3.0* compiled with gcc/acml-4.0.1 (no dotblas); no ScaLAPACK used::

       # p - processes, p0 - reference processes, t - time [sec], s - speedup, e - efficiency
       # GPAW version 6.5092: stages: 1 - initialization, 2 - fixdensity, 3 - SCF, 4 - forces, 5 - total
       # p     p/p0   t1      s1      e1    t2      s2      e2    t3      s3      e3    t4      s4      e4    t5      s5      e5
            16   1.00   190.5    16.0  1.00   823.5    16.0  1.00   563.0    16.0  1.00     0.0     0.0  0.00  1577.0    16.0  1.00
            32   2.00   112.5    27.1  0.85   454.5    29.0  0.91   310.0    29.1  0.91     0.0     0.0  0.00   877.0    28.8  0.90
            64   4.00    71.0    42.9  0.67   255.0    51.7  0.81   172.0    52.4  0.82     0.0     0.0  0.00   498.0    50.7  0.79

   - xeon (HP DL160 G6: 64-bit Intel Nehalem Xeon X5570 quad-core 2.93 GHz CPUs) nodes (ethernet):
     performed on EL5 with gcc43/acml-4.3.0/acml-4.3.0, GPAW **0.6.5092**,
     numpy *1.3.0* compiled with gcc/acml-4.0.1 (no dotblas); no ScaLAPACK used::

       # p - processes, p0 - reference processes, t - time [sec], s - speedup, e - efficiency
       # GPAW version 6.5092: stages: 1 - initialization, 2 - fixdensity, 3 - SCF, 4 - forces, 5 - total
       # p     p/p0   t1      s1      e1    t2      s2      e2    t3      s3      e3    t4      s4      e4    t5      s5      e5
            16   1.00   116.0    16.0  1.00   444.0    16.0  1.00   302.0    16.0  1.00     0.0     0.0  0.00   862.0    16.0  1.00
            32   2.00    66.0    28.1  0.88   270.0    26.3  0.82   184.0    26.3  0.82     0.0     0.0  0.00   520.0    26.5  0.83
            64   4.00    48.0    38.7  0.60   159.0    44.7  0.70   109.0    44.3  0.69     0.0     0.0  0.00   316.0    43.6  0.68

   Clearly SCF part scales better than the initialization stage.
   Using of ScaLAPACK does not result in any noticeable improvement:
   even for the fastest 64 cores run on xeon the diagonalization part
   takes only 4% of the total runtime. This is to be expected from
   a rather small hamiltonian matrix size (1056 bands).
   **Note** that runtimes on opteron ethernet (EL5) and infiniband (EL4) nodes
   are not directly comparable due to different operating system,
   gcc, and numpy versions.

 - for a comparison of what to expect on different machines, the following absolute times where obtained with r=[1,1,1] (without ScaLAPACK)

   ===================   ================= ============  ======= ============  ========  ========      
   host                  type              cpu type      MHz     # procs       time [s]  date
   ===================   ================= ============  ======= ============  ========  ========      
   jump.fz-juelich.de    IBM Regatta p690+ Power4+       1700    2             88        23.3.09
   jump.fz-juelich.de    IBM Regatta p690+ Power4+       1700    4             51        23.3.09
   mmos3                 LINUX             Intel Q6600   2394    2             85        23.3.09
   mmos3                 LINUX             Intel Q6600   2394    4             62        23.3.09
   bfg.uni-freiburg.de   LINUX             Xeon 5160     3000    2             156       23.3.09
   bfg.uni-freiburg.de   LINUX             Xeon 5160     3000    4             119       23.3.09
   ===================   ================= ============  ======= ============  ========  ========      

2. Medium size system
+++++++++++++++++++++

The system used in this benchmark is another one of medium size, as for the year 2008,
and consists of a gold cluster interacting with organic groups
(see `<http://www.pnas.org/cgi/content/abstract/0801001105v1>`_) in a box of 32**3 Angstrom**3,
3366 electrons, and 1728 bands, and 240**3 grid points (grid spacing of ~0.13).
LCAO initialization stage is performed, then 3 SCF steps with a constant
potential and 2 full SCF steps.
All the stages are timed separately, due to their different scaling.

**Note** that the size of the system can be changed easily by modifying
just one variable in :svn:`~doc/devel/Au_cluster/Au_cluster.py`::

  r = [1, 1, 1]

Prerequisites
+++++++++++++

This benchmark requires approximately 2 GB of RAM memory per core
and at least 512 cores, up to 4096.
The amount of disk space required is minimal.

The following packages are required (names given for Fedora Core 10 system):

 - python, python-devel
 - numpy
 - python-matplotlib
 - openmpi, openmpi-devel
 - blacs, scalapack
 - bash
 - `campos-gpaw <https://wiki.fysik.dtu.dk/gpaw/install/installationguide.html>`_
 - `campos-ase3 <https://wiki.fysik.dtu.dk/ase/download.html>`_

**Note** that GPAW has to built with ScaLAPACK enabled -
please refer to :ref:`platforms_and_architectures` for hints on
installing GPAW on different platforms.

Results
+++++++

GPAW code is executed in parallel in order to benchmark a number of processes
that ranges from 256,
through integer powers of 2 and up to the total number of CPU 4096 cores.

The number of bands (1728) and cores are chosen to make comparisons
of different band parallelizations (:ref:`band_parallelization`) possible.

**Note**: to achive optimal performance diagonalization steps are performed
on `5x5` blacs grid with block size of `64` specified by adding ``--gpaw=blacs=1 --sl_diagonalize=5,5,64`` options.

**Note** also that a default domain decomposition is appplied, and different
results can be obtained by tuning ``--domain-decomposition`` argument
to your platform (see :ref:`submit_tool_on_niflheim`).

**Note**: the ``--gpaw=usenewlfc=1`` option is required to skip the calculation of forces
and decrease **memory** usage.

The results of the benchmark is scaling of execution time of different stages
of GPAW run with the number of processes (CPU cores).


Getting the results
+++++++++++++++++++

Please perform the following steps:

 - use the following commands to setup the benchmark::

    bash
    mkdir Au_cluster; cd Au_cluster
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/Au_cluster/Au102_revised.xyz
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/Au_cluster/Au_cluster.py
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/Au_cluster/akka.sh
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/Au_cluster/intrepid.sh
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/Au_cluster/prepare.sh
    wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/256H2O/scaling.py
    # set the prefix directory: results will be in $PATTERN_*_
    export PATTERN=Au_cluster_240_05x05m64.grid
    sh prepare.sh

   **Warning**: the choice of the directory names is not free in the sense that
   the number of processes has to come at the end of directory name,
   and be delimited by two underscores.

 - run with, for example:

    - on akka::

       cd ${PATTERN}_00256_; qsub -l nodes=32:8 ../akka.sh; cd ..
       cd ${PATTERN}_00512_; qsub -l nodes=64:8 ../akka.sh; cd ..
       cd ${PATTERN}_01024_; qsub -l nodes=128:8 ../akka.sh; cd ..
       cd ${PATTERN}_02048_; qsub -l nodes=256:8 ../akka.sh; cd ..
       cd ${PATTERN}_04096_; qsub -l nodes=512:8 ../akka.sh; cd ..

   **Warning**: on Linux clusters it s desirable to repeat these runs 2-3 times
   to make sure that they give reproducible time.

 - analyse the results::

    python scaling.py -v --dir=. --pattern="Au_cluster_240_05x05m64.grid_*_" Au_cluster

   A typical output may look like
   (example given for Intel Xeon dual-socket, quad-core L5k CPUs, 2.5 GHz,
   GPAW linked with Intel mkl, infiniband)::

    # p - processes, p0 - reference processes, t - time [sec], s - speedup, e - efficiency
    # GPAW version 2843: stages: 1 - initialization, 2 - fixdensity, 3 - SCF, 4 - forces, 5 - total
    # p     p/p0   t1      s1      e1    t2      s2      e2    t3      s3      e3    t4      s4      e4    t5      s5      e5
        512   1.00   243.5   512.0  1.00   856.5   512.0  1.00   900.0   512.0  1.00     0.0     0.0  0.00  2000.0   512.0  1.00
       1024   2.00   155.5   801.7  0.78   466.5   940.0  0.92   489.0   942.3  0.92     0.0     0.0  0.00  1111.0   921.7  0.90
       2048   4.00   148.5   839.5  0.41   241.5  1815.9  0.89   256.0  1800.0  0.88     0.0     0.0  0.00   646.0  1585.1  0.77
