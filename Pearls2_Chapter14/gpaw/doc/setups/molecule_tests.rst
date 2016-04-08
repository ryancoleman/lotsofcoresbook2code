.. _molecule_tests:

==============
Molecule tests
==============


**Warning**: this page is outdated.
For performance of GPAW for molecular systems refer to :ref:`g2_1`.

Atomization energies and bond lengths for a set of small molecules
have been calculated with the PBE functional.  All calculations are
done with a grid spacing of 0.16 Å, zero-boundary conditions and
approximately 6 Å of vacuum surrounding the molecules.  Compensation
charges are expanded with correct multipole moments up to
:math:`\ell_{max}=2`.  Open-shell atoms are treated as non-spherical with
integer occupation numbers, and zero-point energy is not included in
the atomization energies. The numbers are compared to very accurate,
state-of-the-art, PBE calculations [1]_.  The script that does the
calculations is :svn:`~gpaw/testing/molecule_test.py`.


Bond lengths
============

.. figure:: bondlengths.png

Bondlengths in Å:

.. csv-table::
   :file: bondlengths.csv		
   :header: dimer, GPAW, reference, error

Atomization energies
====================

Atomization energies in eV:

.. csv-table::
   :file: atomization_energies.csv
   :header: molecule, GPAW, reference, error

References
==========

.. [1] "The Perdew-Burke-Ernzerhof exchange-correlation functional
       applied to the G2-1 test set using a plane-wave basis set",
       J. Paier, R. Hirschl, M. Marsman and G. Kresse,
       J. Chem. Phys. 122, 234102 (2005)

.. [2] "Molecular and Solid State Tests of Density Functional
       Approximations: LSD, GGAs, and Meta-GGAs", S. Kurth,
       J. P. Perdew and P. Blaha, Int. J. Quant. Chem. 75, 889-909
       (1999)

.. [3] "Comment on 'Generalized Gradient Approximation Made Simple'",
       Y. Zhang and W. Yang, Phys. Rev. Lett.

.. [4] Reply to [3]_, J. P. Perdew, K. Burke and M. Ernzerhof
