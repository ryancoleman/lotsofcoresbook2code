.. _jellium:
    
=======
Jellium
=======


In this tutorial, we try to reproduce some old jellium calculations by
Lang and Kohn [Lang70]_.


Bulk
====

Let's do a calculation for `r_s=5` Bohr.  We use a cubic cell with a
lattice constant of `a=1.6` Å, 8*8*8 grid points and and a k-point
sampling of 12*12*12 points.

.. literalinclude:: bulk.py

In the text output from the calculation, one can see that the cell
contains 0.0527907 electrons.


Surfaces
========

Now, we will do a surface calculation.  We put a slab of thickness
16.0 Å inside a box, periodically repeated in the x- and y-directions
only, of size 1.6*1.6*25.6 Å and use 12*12 k-points to sample the
surface BZ:

.. literalinclude:: surface.py

The surface energy is:

>>> from ase.io import read
>>> e0 = read('bulk.txt').get_potential_energy()
>>> e = read('surface.txt').get_potential_energy()
>>> a = 1.6
>>> L = 10 * a
>>> sigma = (e - L / a * e0) / 2 / a**2
>>> print('%.2f mev/Ang^2' % (1000 * sigma))
5.40 mev/Ang^2
>>> print('%.1f erg/cm^2' % (sigma / 6.24150974e-5))
86.6 erg/cm^2

which is reasonably close to the value from Lang and Kohn: 100
ergs/cm\ `^2`.

Here is the electron density profile:

.. literalinclude:: fig2.py

Compare with Fig. 2 in [Lang70]_:

.. image:: fig2.png


Other jellium geometries
========================

For other geometries, one will have to subclass
:class:`~gpaw.jellium.JelliumPoissonSolver`, and implement the
:meth:`~gpaw.jellium.JelliumPoissonSolver.get_mask` method:

.. autoclass:: gpaw.jellium.JelliumPoissonSolver
   :members:


-------------

.. [Lang70] N. D. Lang and W. Kohn,
   Phys. Rev. B 1, 4555-4568 (1970),
   *Theory of Metal Surfaces: Charge Density and Surface Energy*,
   DOI: 10.1103/PhysRevB.1.4555

