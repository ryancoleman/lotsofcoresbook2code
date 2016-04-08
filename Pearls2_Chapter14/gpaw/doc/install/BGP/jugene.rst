.. _jugene:

====================
jugene.fz-juelich.de
====================

Here you find information about the system
`<http://www.fz-juelich.de/jsc/jugene>`_.

Numpy needs to be build with powerpc-bgp-linux-gfortran instead of gfortran
compiler, so in order to build numpy specify the environment variable F90::

 $ export F90=/bgsys/drivers/ppcfloor/gnu-linux/powerpc-bgp-linux/bin/gfortran

After that, numpy can be installed to $HOME/python as::

 $ ldpath=/bgsys/drivers/ppcfloor/gnu-linux/lib
 $ p=/bgsys/drivers/ppcfloor/gnu-linux/bin/python
 $ LD_LIBRARY_PATH="$ldpath" $p setup.py install --home=$HOME/python

In order to build GPAW, use the following customize.py::

 scalapack = True

 libraries = [
            'scalapack',
            'blacsCinit',
            'blacsF77init',
            'blacs',
            'lapack',
            'esslbg',
             'xl',
             'xlopt',
             'xlf90_r',
             'xlfmath',
             'pthread',
             'xlomp_ser',
            ]

 library_dirs = [
           '/bgsys/local/lib/',
            '/opt/ibmcmp/xlf/bg/11.1/lib',
            '/opt/ibmcmp/xlsmp/bg/1.7/lib',
            '/bgsys/drivers/ppcfloor/gnu-linux/lib'
            ]

 extra_compile_args += ['-std=c99']

 define_macros += [('GPAW_AIX', '1')]
 define_macros += [('GPAW_BGP', '1')]
 define_macros += [('GPAW_NO_UNDERSCORE_BLAS', '1')]
 define_macros += [('GPAW_NO_UNDERSCORE_LAPACK', '1')]
 define_macros += [('GPAW_NO_UNDERSCORE_CBLACS', '1')]
 define_macros += [('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]
 define_macros += [('GPAW_NO_UNDERSCORE_BLACS', '1')]
 define_macros += [('GPAW_NO_UNDERSCORE_SCALAPACK', '1')]

Because of missing ``popen3`` function you need to remove all the
contents of the :file:`gpaw/version.py` file after ``version =
'0.4'``.  The same holds for :file:`ase/version.py` in the ase
installation!  Suggestions how to skip the ``popen3`` testing in
:file:`gpaw/version.py` on BGP are welcome!


Here is an example of batch job script::

  #!/bin/bash
 # @ job_name = hello
 # @ output = $(job_name).o$(jobid)
 # @ error = $(job_name).e$(jobid)
 # @ wall_clock_limit = 00:12:00
 # @ notification = never
 # @ notify_user = my_email@csc.fi
 # @ job_type = bluegene
 # @ bg_size = 1
 # @ queue
 home=/homea/prace/prace025
 prog=${home}/python/bin/gpaw-python
 args="${home}/test-gpaw/CH4.py"

 mpirun=/bgsys/drivers/ppcfloor/bin/mpirun

 ldpath="/bgsys/local/lib/ibmcmp/lib/bglib"

 pythonpath="${home}/python/lib/python/"

 gpaw_setups="${home}/gpaw-setups-0.4.2039"

 runargs="-np 4"

 runargs="$runargs -cwd $PWD"
 runargs="$runargs -mode SMP"
 runargs="$runargs -env LD_LIBRARY_PATH=$ldpath -env PYTHONPATH=$pythonpath -envGPAW_SETUP_PATH=$gpaw_setups"

 echo "Hello. This is `hostname` at `date` `pwd`"

 echo "$mpirun $runargs -exe $prog $args"
 /usr/bin/time $mpirun $runargs -exe $prog -args $args

 echo "Program completed at `date` with exit code $?."

The batch jobs are submitted with ``llsubmit``::

 $ llsubmit job_file
