======================
Adding new calculators
======================

Adding an ASE interface to your favorite force-calculator is very
simple.  Take a look at the
:class:`~ase.calculators.calculator.Calculator` and
:class:`~ase.calculators.calculator.FileIOCalculator` classes below
(the code is here: :trac:`ase/calculators/calculator.py`).  You should
inherit from the :class:`~ase.calculators.calculator.FileIOCalculator`
and implement the :meth:`~ase.calculators.calculator.Calculator.read`,
:meth:`~ase.calculators.calculator.FileIOCalculator.read_results` and
:meth:`~ase.calculators.calculator.FileIOCalculator.write_input` methods.
The methods :meth:`~ase.calculators.calculator.Calculator.set`,
:meth:`~ase.calculators.calculator.Calculator.check_state` and
:meth:`~ase.calculators.calculator.Calculator.set_label` may also need
to be implemented.

.. seealso::

   * The code for our Abinit interface: :trac:`ase/calculators/abinit.py`
   * :ref:`aep1`
   * :mod:`ase.calculators`


Description of base-classes
===========================


The Calculator base-class
-------------------------

.. autoclass:: ase.calculators.calculator.Calculator
   :members:
   :private-members:
   :member-order: bysource


The FileIOCalculator class
--------------------------

.. autoclass:: ase.calculators.calculator.FileIOCalculator
   :members:
   :private-members:
   :member-order: bysource
