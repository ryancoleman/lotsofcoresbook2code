.. _mira:

=======================
Blue Gene/Q - Mira
=======================

The build instructions here are representative of the Blue Gene/Q
system at the Argonne Leadership Computing Facility. Users will need
to adapt these instructions for specific installations at their
respective sites.

In addition to the standard libraries needed on other platforms,
Scalable Python `<https://gitorious.org/scalable-python>`_ is required
for running effectively on Blue Gene/Q. A build script for Scalable
Python is provided below:

.. literalinclude:: build_scalable_python.sh

NumPy 1.3.0 or later is recommended. Disutils does not work well on
PPC architectures and a compiler must be explictly specified. A build
script for NumPy 1.3.0 is provided below:

.. literalinclude:: build_numpy.sh

GPAW will build with the XL legacy MPI wrapper script. It is
recommeded that you statically link as many libraries as possible into
GPAW to avoid potential performance bottlencks in loading shared
libraries at scale. This can be done with some modification of the
stock GPAW config.py file :svn:`~doc/install/BGQ/config_mira.py`.

Lastly, we recommend that GPAW is compiled with both ScaLAPACK
(v. 2.0.2 or later) as well as HDF5 support. Here is an example
customization file:

.. literalinclude:: customize_mira_xlc_mpi.py

which requires a number of wrappers for the XL compilers
(:svn:`~doc/install/BGQ/bgq_xlc.py` and
:svn:`~doc/install/BGQ/bgq_xlc_linker.py`). A build script for GPAW
is provided for convenience :svn:`~doc/install/BGQ/build_gpaw.sh`.

After all Python modules are installed, they should be byte compiled
before running GPAW. This can be accomplished by going to the top level
directory for each Python library (Scalable Python, NumPy, ASE, and
GPAW) an executing the command::

  ${python} -m compileall .

where ``${python}`` is the explicit path to the Scalable Python
interpreter.

Some users have noticed that the Python interpreter may unnecessarily
re-compile Python modules. This is problematic at large (>10,000)
number of MPI tasks and we recommend that users set the environment
variable::

  PYTHONDONTBYTECOMPILE=1

in the job submissions script.









