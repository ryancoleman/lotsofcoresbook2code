.. _structure_optimization:

========================
 Structure optimization
========================

In the tutorial on :ref:`how to calculate atomization energies
<atomization_energy>`, we calculated the atomization energy for
:math:`\rm{H}_2` using the experimental bond length of 0.74 Å.  In
this tutorial, we ask a :ase:`QuasiNewton
<ase/optimize.html#module-optimize.qn>` minimizer to iteratively find
the structural energy minimum, where all atomic forces are below 0.05
eV/Å.  The following script will do the job:

.. literalinclude:: relax.py

The result is:

.. literalinclude:: optimization.txt

.. note::
   You must run the :ref:`atomization <atomization_energy>` script first.

To save time you could have told the minimizer to keep one atom fixed,
and only relaxing the other. This is achieved through the use of
constraints::

  molecule.set_constraint(FixAtoms(mask=[0, 1]))

The keyword ``mask`` contains list of booleans for each atom indicating
whether the atom's position should be fixed or not. See the
:ase:`constraints <ase/constraints.html>` section on the ASE page for
more information and examples for setting constraints.
