==========================
Dipole corrections in GPAW
==========================

As an example system, a 2 layer :math:`2\times2` slab of fcc (100) Al
is constructed with a single Na adsorbed on one side of the surface.

.. literalinclude:: dipole.py
    :lines: 1-14

.. image:: slab.png

The :func:`ase.lattice.surface.fcc100` function will create a slab
with periodic boundary conditions in the xy-plane only and GPAW will
therefore use zeros boundary conditions for the the wave functions and
the electrostatic potential in the z-direction as shown here:

.. image:: zero.png

The blue line is the xy-averaged potential and the green line is the
fermi-level.

.. note::

    You need a bit of magic to get the electrostatic potential from a
    gpw file:

    >>> from ase.units import Hartree
    >>> from gpaw import GPAW
    >>> calc = GPAW('zero.gpw', txt=None)
    >>> calc.restore_state()
    >>> v = calc.hamiltonian.vHt_g * Hartree
    >>> v.shape
    (56, 56, 167)

If we use periodic boundary conditions in all directions:

.. literalinclude:: dipole.py
    :lines: 16-19

the electrostatic potential will be periodic and average to zero:

.. image:: periodic.png

In order to estimate the work functions on the two sides of the slab,
we need to have flat potentials or zero electric field in the vacuum
region away from the slab.  This can be achieved by using a dipole
correction:

.. literalinclude:: dipole.py
    :lines: 21-25

.. image:: corrected.png

.. warning::

    * Information about use of a dipole correction is currently not
      written to the gpw file.  See below how to restart such a
      calculation.

See the full Python script here: :download:`dipole.py`.  The script
used to create the figures in this tutorial is shown here:

.. literalinclude:: plot.py

.. autoclass:: gpaw.dipole_correction.DipoleCorrection
    :members:
