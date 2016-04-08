.. module:: ase.utils

==============================
Utillity functions and classes
==============================

.. autofunction:: ase.utils.opencew
.. autofunction:: ase.utils.gcd

Timing is facilitated with the :class:`~ase.utils.timing.Timer` class
that can be accessed most easily by the 
:class:`~ase.utils.timing.timer` class:

.. autoclass:: ase.utils.timing.timer

.. index:: Bulk modulus

Equation of state
=================

The :class:`~ase.utils.eos.EquationOfState` class can be used to find
equilibrium volume, energy, and bulk modulus for solids:

.. autoclass:: ase.utils.eos.EquationOfState
  :members: fit, plot


.. seealso::  The :ref:`eos` tutorial.


Symmetry analysis
=================

http://spglib.sourceforge.net/pyspglibForASE/


Phonons
=======

http://phonopy.sourceforge.net/
