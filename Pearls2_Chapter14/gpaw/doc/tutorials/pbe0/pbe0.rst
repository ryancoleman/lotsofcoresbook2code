.. _pbe0_tut:

==================================
PBE0 calculations for bulk silicon
==================================

.. seealso::

   * :ref:`bandstructures` tutorial.
   * :ref:`band exercise` exercice.


PBE band structure
==================

GPAW can currently not do self-consistent PBE0 calculations, so this
tutorial will do non-selfconsistent PBE based on selfconsistent PBE.
Here is a little function that will calculate the PBE groundstate:

.. literalinclude:: si_pbe.py

We do a calculation for a lattice constant of 5.43 Å and a k-point
sampling of 8*8*8 points::

    si = groundstate(5.43, 8)

and write the gpw-file *including* wave functions for later use::

    si.calc.write('Si-PBE.gpw', mode='all')


.. tip::

    Alternative for command-line folks:  Use the :program:`gpaw` command.

    .. highlight:: bash

    ::
      
        $ gpaw Si -x diamond -a 5.43 --kpts=8,8,8 --xc=PBE --write-gpw-file=all -i PBE

    .. highlight:: python


The band structure can be calculated like this:

.. literalinclude:: bs.py

and plotted like this:

.. literalinclude:: bs_plot.py

.. image:: bs-PBE.png


PBE0 band structure
===================

In order to plot a nice band structure for PBE0, we are forced to
interpolate between the 8*8*8 k-points.  We do that with Wannier
functions.  The initial guess for 4 localized Wannier will be 4
spherical gaussians located on the 4 Si-Si bonds in the unit cell:

.. literalinclude:: wannier.py

Read in the rotation matrices and perform the interploation for both
the PBE and the non self-consistent PBE0 eigenvalues:

.. literalinclude:: bs_pbe0.py

and plot:

.. literalinclude:: bs_plot2.py

.. image:: bs-PBE0.png


Special points
--------------

TODO XXX ...


Lattice constant and bulk modulus
=================================

.. image:: a.png

.. image:: B.png
