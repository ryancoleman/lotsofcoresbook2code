.. _neb:

=========================================
NEB calculations parallelized over images
=========================================

The :ase:`Gold atom diffusion on Al(100)
<tutorials/neb/diffusion.html>` example can be used with GPAW like
this:

.. literalinclude:: neb.py

If we run the job on 12 cpu's::

  $ mpirun -np 12 gpaw-python neb.py

then each of the three internal images will be parallelized over 4 cpu's.

The results are read with::

  $ ase-gui -n -1 neb?.traj

The energy barrier is found to be ~0.14 eV.
