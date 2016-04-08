.. _building_with_xlc_on_surveyor:

==========================
Building with xlc compiler
==========================

NumPy
=======

We currently do not know how to build NumPy with xlc on BG/P.

GPAW
====

A performance improvement of 25% has been observed using xlc over gcc
for medium-size systems where DGEMM does not dominate the wall-clock
time. For large-systems, NBANDS > 5000, there is is little performance 
improvement since DGEMM makes up a large fraction of the wall-clock
time.

Proceed as in the :ref:`building_with_gcc_on_surveyor`,
but use the following :svn:`bgp_xlc.py` file:

.. literalinclude:: bgp_xlc.py

Finally, change the lines in :svn:`customize_surveyor_gcc.py` accordingly::

  mpicompiler = "bgp_xlc.py"
  compiler = "bgp_xlc.py"
  mpilinker = "bgp_xlc_linker.py"

Everything else should be the same.

