===============================
Setting up an aluminium surface
===============================

In this exercise, we build an Al(100) surface. For this
surface, we calculate the surface energy and other properties.



Fcc Surfaces
============

A real fcc surface has a large number of atomic layers, but most
surface properties are well reproduced by a slab with just 2 - 20
layers, depending on the material and what properties you are looking
for.

The most important cubic surfaces are (100), (110), and (111).  For
face centered cubic, (111) has the most compact atomic arrangement,
whereas (110) is most open. Here we'll focus on (100).

* What is the coordination number *Z* (number of nearest neighbors) of
  an fcc(100) surface atom?  What is it for a bulk atom?  Start the
  Python interpreter and try this::

    from ase.visualize import view
    from ase.lattice.surface import fcc100
    s = fcc100('Al', (1, 1, 5))
    view(s, repeat=(4, 4, 1))

  Read more about the :func:`~ase.lattice.surface.fcc100` function
  :ase:`here <ase/lattice.html>`.

* Answer the same questions for the (111) surface.


Aluminum fcc(100) surface energetics
====================================

One surface property is the surface tension `\sigma` defined
implicitly via:

.. math:: E_N = 2A\sigma + NE_B

where `E_N` is the total energy of a slab with `N` layers,
`A` the area of the surface unit cell (the factor 2 because the slab
has two surfaces), and finally `E_B` is the total energy per bulk
atom.  The limit `N \rightarrow \infty` corresponds to the macroscopic
crystal termination.

Estimate the surface tension using an expression from the simplest
Effective Medium Theory (EMT) description:

.. math:: A\sigma \simeq [1 - (Z/Z_0)^{1/2}] E_{coh}

where `Z` and `Z_0` are the coordination numbers (number of nearest
neighbors) of a surface and a bulk atom, respectively, and `A` is the
surface area per surface atom, and `E_{coh} = E_{atom}-E_B > 0` is
the cohesive energy per bulk atom. For Aluminium we have `E_{coh}` = 3.34 eV.

* Derive the following equation:

  .. math:: \sigma = \frac{NE_{N-1} - (N-1)E_N}{2A}

* The script :svn:`~doc/exercises/surface/Al100.py` defines
  an ``energy()`` function for calculating `E_N`.  Use it to calculate
  `\sigma` for `N` = 5.  Use a two-dimensional Monkhorst-Pack
  **k**-point sampling (``kpts=(k, k, 1)``) that matches the size of
  your unit cell.  

  .. hint::

    A rule of thumb for choosing the initial **k**-point sampling is,
    that the product, *ka*, between the number of **k**-points, *k*,
    in any direction, and the length of the basis vector in this
    direction, *a*, should be:

    * *ka* ~ 30 Å, for *d* band metals
    * *ka* ~ 25 Å, for simple metals
    * *ka* ~ 20 Å, for semiconductors
    * *ka* ~ 15 Å, for insulators

    Remember that convergence in this parameter should always be checked.

* How well is the EMT estimate satisfied?

.. note:: The experimental value of `\sigma` is 54 meV/Å\ :sup:`2`. 
   However, this was obtained from the curvature of an aluminium drop and
   is more likely to represent the value for a closepacked Al(111) surface.


Work function
=============

Run the :svn:`~doc/exercises/surface/work_function.py`
script and estimate the work function for a Al(100) surface (this
script does not run in parallel). A typical
experimental value for the work function of the Al(100) surface is
4.20 eV. 
