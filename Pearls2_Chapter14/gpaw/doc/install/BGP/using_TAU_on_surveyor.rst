.. _using_TAU_on_surveyor:

=====================
Using TAU on surveyor
=====================

Start by reading the `main profiling page <https://wiki.fysik.dtu.dk/gpaw/devel/profiling.html>`_

Do **not** using the customize.py from the above page, following the instructions found on this
page instead. The following mostly applies to automatic instrumentation using the TAU compiler
scripts. 

Overhead with the TAU for automatic instrumentation has been measured at about 20% for the
`b256H2O.py <https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/256H2O/b256H2O.py>`_.
Use manual instrumentation if this overhead is unacceptable or you might need to 
add more functions to the existing selective instrumentation file
`select.tau <https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/profiling/select.tau>`_

The version of TAU at ALCF is updated frequently, please check:
`<https://wiki.alcf.anl.gov/index.php/Tuning_and_Analysis_Utilities_(TAU)>`_

We first described the case of automatically instrumentation with the
GCC and XLC compilers.

Automatic Instrumentation
============================
The instructions below are for using source-based automatic instrumentation with TAU. It is well tested with the GCC compiler. It should also work with XLC compiler, but is untested.

Add the following lines to you ``.softenvrc``::

  TAUARCHITECTURE = bgp
  TAUVERSION = 2.19.2
  TAU_MAKEFILE = /soft/apps/tau/tau-$TAUVERSION/bgp/lib/Makefile.tau-bgptimers-gnu-mpi-python-pdt
  TAU_OPTIONS = '-optVerbose -optShared -optTauSelectFile=select.tau \
  	      -optTau="-rn Py_RETURN_NONE -i/soft/apps/tau/tau-'$TAUVERSION'/include/TAU_PYTHON_FIX.h"'
  PATH += /soft/apps/tau/tau-$TAUVERSION/$TAUARCHITECTURE/bin

The bindings are located in
``/soft/apps/tau/tau/tau-$TAUVERSION/bgp/lib/<name>``.  This particular TAU library binding supports BGP timers (a low-level
timer with minimal overhead), MPI, GNU compiler, and Python. This is the recommended library binding for flat profiles.

You will need to change one line in your :svn:`bgp_gcc.py`::

  cmd = "tau_cc.sh %s %s"%(flags, cmd)
  
It is recommend that dynamic linking be used for the MPI libraries
instead of static linking, see the comments in :svn:`configure_surveyor.py`

Additionally, the TAU library must be linked manually. See the
comments in :svn:`customize_surveyor.py`

Manual Instrumentation
=============================
It is straightforward to swap the default Python timers with those
provided by TAU and make subclass of the default GPAW calculator::

  from gpaw.utilities.timing import TAUTimer
 
  class MyGPAW(GPAW):
         timer_class = TAUTimer

  calc = MyGPAW(...)

Run time environment variables
================================
Please see:
`<https://wiki.alcf.anl.gov/index.php/Tuning_and_Analysis_Utilities_(TAU)#Running_With_TAU>`_.

Here are the recommended run time environment variables that should be passed to Cobalt via qsub::

  TAU_VERBOSE=1:TAU_THROTTLE=0:TAU_COMPENSATE=1:TAU_METRICS=BGPTIMERS

TAU_COMPENSATE seems to cause problems with manual instrumentation, so do not set it to 0 which
means off. In any case, it should not be particularly relevant unless you have manual timers on a
frequently accessed lightweight functions.

Submitting jobs
==================

Modification to the submission script for use with TAU are show in the
example submission script :svn:`surveyor.sh`

If you are doing manual instrumentation, simply pass the actual input file to ``gpaw-python`` instead. For automatic instrumentation, you need to ``wrapper.py`` instead::

  import tau

  def OurMain():
      import CH4;

  tau.run('OurMain()')

TAU run will then produce ``profile.*`` files that can be merged into
the default TAU's ``ppk`` format using the command issued from the directory
where the ``profile.*`` files reside::

 paraprof --pack CH4.ppk

The actual analysis can be made on a different machine, by transferring
the ``CH4.ppk`` file from ``surveyor``, installing TAU, and launching::

 paraprof CH4.ppk
