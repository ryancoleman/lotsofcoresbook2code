.. _diffusion_tutorial:

===============================================
Diffusion of gold atom on Al(100) surface (NEB)
===============================================

First, set up the initial and final states:

|initial| |final|

.. literalinclude:: diffusion1.py

.. note::  Notice how the tags are used to select the constrained atoms

Now, do the NEB calculation:

.. literalinclude:: diffusion2.py

Visualize the results with::

   ase-gui neb.traj 

and select Tools->NEB.

|ts| |barrier|

.. note::

   For this reaction, the reaction coordinate is very simple: The
   *x*-coordinate of the Au atom.  In such cases, the NEB method is
   overkill, and a simple constraint method should be used like in this
   tutorial: :ref:`constraints_diffusion_tutorial`.

.. seealso::

   * :mod:`ase.neb`
   * :mod:`ase.constraints`
   * :ref:`constraints_diffusion_tutorial`
   * :func:`~ase.lattice.surface.fcc100`
   


.. |initial| image:: diffusion-I.png
.. |final| image:: diffusion-F.png
.. |ts| image:: diffusion-T.png
.. |barrier| image:: diffusion-barrier.png


Restarting NEB
==============

Restart NEB from the trajectory file:

.. literalinclude:: diffusion4.py


Parallelizing over images with MPI
==================================

Instead of having one process do the calculations for all three
internal images in turn, it will be faster to have three processes do
one image each. In order to be able to run python with MPI
you need a special parallel python interpreter, for example gpaw-python.

The example below can then be run
with ``mpiexec -np 3 gpaw-python diffusion3.py``:

.. literalinclude:: diffusion3.py


