.. _profiling:

=========
Profiling
=========

profile
=======

Python has a ``profile`` module to help you find the places in the
code where the time is spent.

Let's say you have a script
(`CH4.py <https://svn.fysik.dtu.dk/projects/gpaw/trunk/test/CH4.py>`_)
that you want to run through the profiler.  This is what you do:

>>> import profile
>>> profile.run('import CH4', 'prof')

This will run your script and generate a profile in the file ``prof``.
You can also generate the profile by inserting a line like this in
your script::

  ...
  import profile
  profile.run('atoms.get_potential_energy()', 'prof')
  ...

To analyse the results, you do this::

 >>> import pstats
 >>> pstats.Stats('prof').strip_dirs().sort_stats('time').print_stats(20)
 Tue Oct 14 19:08:54 2008    prof

         1093215 function calls (1091618 primitive calls) in 37.430 CPU seconds

   Ordered by: internal time
   List reduced from 1318 to 20 due to restriction <20>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
    37074   10.310    0.000   10.310    0.000 :0(calculate_spinpaired)
     1659    4.780    0.003    4.780    0.003 :0(relax)
   167331    3.990    0.000    3.990    0.000 :0(dot)
     7559    3.440    0.000    3.440    0.000 :0(apply)
      370    2.730    0.007   17.090    0.046 xc_correction.py:130(calculate_energy_and_derivatives)
    37000    0.780    0.000    9.650    0.000 xc_functional.py:657(get_energy_and_potential_spinpaired)
    37074    0.720    0.000   12.990    0.000 xc_functional.py:346(calculate_spinpaired)
      ...
      ...

The list shows the 20 functions where the most time is spent.  Check
the pstats_ documentation if you want to do more fancy things.

.. _pstats: http://docs.python.org/lib/module-profile.html


.. tip::

   Since the ``profile`` module does not time calls to C-code, it
   is a good idea to run the code in debug mode - this will wrap
   calls to C-code in Python functions::

     $ python script.py --debug

.. tip::

   There is also a quick and simple way to profile a script::

     $ pyhton script.py --profile=prof

   This will produce a file called ``prof.0000`` where ``0000`` is the
   rank number (if your run the script in parallel, there will be one
   file per rank).

   Use::

     $ python script.py --profile=-

   to write the report directly to standard output.




TAU
===

`TAU Performance System <http://www.cs.uoregon.edu/research/tau/>`_
is a portable profiling and tracing toolkit for performance analysis
of parallel programs written in Fortran, C, C++, Java, and Python.

Installation and configuration
------------------------------

TAU requires `Program Database Toolkit
<http://www.cs.uoregon.edu/research/pdt/>`_ to produce profiling
data. Follow the vendor installation instructions for ``pdtoolkit``.
**Note**: if you want to use only TAU's ``paraprof`` and/or ``perfexplorer`` for
analysing profile data made elsewhere - skip the installation of
``pdtoolkit``.
Follow the vendor installation instructions for ``TAU``.

PerfDMF are a set of utilities that allow you to create a database for storing and analayzing your profile data. If you plan to collect any significant amount of profile data, it is highly recommend that you set one up. There are a number of options for the perfmdf database. For must users, the default based on derby is the simplest and it can be created by ``perfdmf_configure --create-default``. Other options are
also available:
`<http://www.cs.uoregon.edu/research/tau/docs/newguide/ch18s02.html>`_
Note that altough Paraprof is the simplest way to load profile data into your database,
it is also possible to accomplish this directly with the PerfDMF utilities.
`<http://www.cs.uoregon.edu/research/tau/docs/newguide/ch19.html>`_

Configure the PerfDMF environment by running first ``perfdmf_configure``. 
Choose the default answer for most of the questions, the only important 
points are to provide a reasonable path to the database directory 
(the database may grow to GB's), and ignore the password question:: 
	 
  Please enter the path to the database directory. 
  (/home/camp/user/.ParaProf/perfdmf):/scratch/user/perfdmf 
  Please enter the database username. 
  (): 
  Store the database password in CLEAR TEXT in your configuration file? (y/n):y 
  Please enter the database password: 

  Would you like to upload the schema? [y/n]: y
 
Similarly, run ``perfexplorer_configure`` letting the default settings::

  Would you like to attempt to automatically download the Weka jar file? (y/n) y

  Would you like to attempt to automatically download the required jar files? (y/n) y

  Perfexplorer tables not found.Would you like to upload the schema? [y/n]: y

Generating profile data
------------------------
TAU has a number of capabilities including generating a call path, memory profiling, measuring MPI message sizes, and much more. Here we describe the flat profile. It is the most basic type of profile. It will show you where the wall-clock time is going. There are two methods for generating a flat profile:

Manual
^^^^^^^^^^
This is the most intuitive way to create profiles however it requires that you understand gpaw's program flow. A TauTimer class is available that includes the most essential profiling commands. See gpaw/utilities/timing.py. Time is measured  **only** for each instance of ``timer.start(<text>)`` and ``timer.stop(<text>)``. There are a number of pre-defined timers for the most time consuming parts of GPAW, e.g. RMM-DIIS, subspace diagonalization, etc. It is very straightforward to add your own timers.

Here is an example of how the TauTimer class can be used to profile a calculation. Note that the TAU binding library (\*.so) **must** be in your *PYTHONPATH* and *LD_LIBRARY_PATH*::

  from gpaw.utilties.timing import TauTimer

  class MyGPAW(GPAW):
         timer_class = TauTimer

  calc = MyGPAW(<args>)


Automatic
^^^^^^^^^^^^
Timing information for every Python and C function is measured. You will need to compile a special version of gpaw. This is often referred to as the instrumented binary. It is important to understand that instrumentation will be performed on three distinct levels:
* MPI
* Python
* C

**Note**: instructions were tested with pdtoolkit-3.14.1 and tau-2.18.2p4.
It is necessary to remove any mpi libraries from the *libraries* variable
in ``customize.py``: TAU will perform linking on it's own.

Simply add the following to the ``customize.py`` and run ``python setup.py build_ext --remove-default-flags``::

  import tau
  import os
  tau_path = tau.__file__[0:tau.__file__.find('lib')]
  tau_make = tau_path+'lib/Makefile.tau-mpi-pthread-python-pdt'
  mpicompiler = "tau_cc.sh -tau_options='-optShared -optCompInst -optVerbose -optMpi' -optTau='-rn Py_RETURN_NONE -i"+os.path.join(os.environ['TAUROOT'], 'include', 'TAU_PYTHON_FIX.h')+"' -tau_makefile="+tau_make
  #mpicompiler = "tau_cc.sh -tau_options='-optShared -optCompInst -optVerbose -optMpi' -optTau='-rn Py_RETURN_NONE' -tau_makefile="+tau_make
  mpilinker = mpicompiler
  compiler = mpicompiler

  extra_link_args += ['-Wl,-rpath='+tau_path+'lib/']

There may be a number of Makefile TAU stubs available. Choose the one that is appropriate for the profile data that you wish to collect and the compiler. Because automatic instrumentation generally has larger overhead than manual instrumentation, it is recommended to set the compensate option. In this way, the instrumentation time will be substracted out from the time reported by TAU. Without this compensation option, light weight functions may be over-represented in the flat profile.

You may set the following::

  export TAU_VERBOSE=1
  export TAU_THROTTLE=0
  export TAU_COMPENSATE=1

Note that an alternate (and simpler) way to specify the long parameter list to ``tau_cc.sh`` is through the use of the environment variables *TAU_MAKEFILE* and *TAU_OPTIONS*. Ultimately, the the profile data collected and how it is collected is determined by these two environment variables.
 
To obtain the profiler data run the following ``wrapper.py``::

  import tau

  def OurMain():
      import CH4;

  tau.run('OurMain()')

e.g., for two processes::

  mpirun -np 2 gpaw-python wrapper.py

This will generate ``profile.?.?.?`` files, convert
these files into a ppk (ParaProf Packed Profile) file with::

  paraprof --pack CH4.ppk

You should be able to quickly view the profiler data with::

  paraprof CH4.ppk

Understanding TAU_OPTIONS
---------------------------
There are a number of other *TAU_OPTIONS* which are helpful but
may not work if TAU is not configured correctly.

* **-optCompInst**: Performs the compiler-based instrumentation. This enables instrumentation by modifying the object files. This is only supported with certain compilers. The default is source-based instrumentation. This applies only to C code.
* **-optShared**: Specifies the use of a dynamic TAU library binding (\*.so) instead of the default static library that would otherwise be linked into ``gpaw-python``. Instrumentation on all three distinct levels requires this option. This should **always** be used, otherwise profile information will only be collected for the Python layer. Additionally,  the library binding is chosen at runtime by specifying the TAU library binding directory in your *PYTHONPATH* and *LD_LIBRARY_PATH*
* **-optTau**: This is frequently very platform and compiler specific. See the BG/P page for examples.
* **-optTauSelectFile**: A file containing a list of functions that are to be excluded or included. This is necessary to reduce the instrumentation overhead from lightweight function calls. The following selective instrumentation file is recommended for use with TAU as a number of functions called by libxc created substantial overhead: `select.tau <https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/profiling/select.tau>`_
* **-optVerbose**: Useful for debugging, all the details of the invocation of ``tau_cc.sh`` are passed to stdout.

Analysing profile data
-----------------------

Now, assuming you have an ppk (ParaProf Packed Profile) file ready,
run ``paraprof`` and choose the following using right clicks:
``Applications -> Default -> Add application -> Add experiment -> Add
trial -> Trial Type: ParaProf Packed Profile``.

``paraprof`` allows you to investigate profiler data for a single run (trial).
Repeat the previous step (adding a trial) for parallel runs
with increasing number of processes, exit ``paraprof`` (derby database
format can be accessed by only one program at a time), and run
``perfexplorer`` to investigate the strong scaling of your application.
