.. _huygens:

===============================================
huygens.sara.nl  (IBM Power6, Infiniband, ESSL)
===============================================

Here you find information about the the system
`<http://huygens.supercomputer.nl/description/>`_.

One should not use the systems defaul python, but load the python module::

 $ module load python

Now, numpy and ASE can be installed in the standard way::

 $ python setup.py install ...

In order to use gcc for parallel compilation of GPAW, set the environment variable 
MP_COMPILER::
 
 $ export MP_COMPILER=gcc

Use the following customize.py::

 libraries += ['xlf90_r', 'xlsmp', 'xlfmath', 'lapack', 'essl', 'xl']

 library_dirs += ['/sara/sw/lapack/3.1.1/lib',
                  '/opt/ibmcmp/xlf/12.1/lib64/',
                  '/opt/ibmcmp/xlsmp/1.8/lib64/',
                 ]

 define_macros += [('GPAW_AIX', '1')]
 define_macros += [('GPAW_NO_UNDERSCORE_BLAS', '1')]
 define_macros += [('GPAW_NO_UNDERSCORE_LAPACK', '1')]
 extra_compile_args += ['-std=c99']
 mpicompiler = 'mpcc'
 mpilinker = 'mpcc

Here is an example batch job script ::

 # @ node = 1
 # @ tasks_per_node = 4
 #
 # Loadleveler can send email, but for this job, we ask Loadleveler not
 #     to send any email:
 #
 # @ notification = never
 #
 # Define the standard input, output and error for the job:
 #
 # @ input = /dev/null
 # @ output = out.$(jobid)
 # @ error = err.$(jobid)
 #
 # @ wall_clock_limit = 0:30:00
 #
 # @ job_type = parallel
 #
 # @ network.MPI = sn_all,not_shared,US
 # @ queue
 #

 cd $HOME/gpaw-benchmarks/

 export PYTHONPATH=$HOME/python/lib/python
 export GPAW_SETUP_PATH=$HOME/gpaw-setups-0.4.2039
 export GPAW_PYTHON=$HOME/python/bin/gpaw-python
 $GPAW_PYTHON input.py


The batch jobs are submitted with ``llsubmit``::

 $ llsubmit job_file

