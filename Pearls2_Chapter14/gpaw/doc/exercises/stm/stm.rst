===============
STM simulations
===============

Scanning Tunneling Microscopy (STM) is a widely used experimental
technique. STM maps out a convolution of the geometric and electronic
structure of a given surface and it is often difficult if not
impossiple to intrepret STM images without the aid of theoretical
tools.

We will use GPAW to simulate an STM image.  Start by doing an Al(100)
surface with hydrogen adsorbed in the ontop site:
:svn:`~doc/exercises/stm/HAl100.py`.  This will produce a
:file:`gpw` file containing the wave functions that are needed for
calculating local density of states.

The STM image can be calculated with the
:svn:`~doc/exercises/stm/stm.py` script::

  $ python stm.py HAl100.gpw

Try the following:

* clean slab without hydrogen
* different number of layers
* different number of **k**-points
* different current
* negative bias (positive tip)
