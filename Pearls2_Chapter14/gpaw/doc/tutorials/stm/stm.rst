.. _stm tutorial:

==============================
Tutorial: STM images - Al(111)
==============================

Let's make a 2 layer Al(111) fcc surface using the
:func:`ase.lattice.surface.fcc111` function:

.. literalinclude:: al111.py
   :lines: 1-3

Now we calculate the wave functions and write them to a file:

.. literalinclude:: al111.py
   :lines: 5-


2-d scans
=========

First initialize the :class:`~ase.dft.stm.STM` object and get the
averaged current at `z=8.0` Å (for our surface, the top layer is at `z=6.338`
Å):

.. literalinclude:: stm.py
   :end-before: matplotlib

From the current we make a scan to get a 2-d array of constant current
height and make a contour plot:

.. literalinclude:: stm.py
   :start-after: scan
   :end-before: figure

.. image:: 2d.png


Linescans
=========

Here is how to make a line-scan:

.. literalinclude:: stm.py
   :start-after: 2d.png

.. image:: line.png
