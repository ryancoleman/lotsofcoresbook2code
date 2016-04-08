.. _scf_conv_eval:

==========================
SCF Convergence Evaluation
==========================
This page presents influence of the choice of different parameters on
the number of steps in the self-consistent total energy calculation.

G2-1
====
PBE0 total energy calculation of molecules from the G2-1 set.
All systems are calculated in spin-polarized mode, with the default
number of bands (left unspecified), unless stated otherwise.

On the plots: the number of self-consistent steps (stats), the total run time,
the systems with the largest number of steps, and the number of
systems for which convergence succeded.

In the corresponding tables: the most common value of the total energy
("energy"), and the differences (run - "energy").
The dot sign denotes the above difference below a threshold
(same as the printed precision of "energy" in eV),
and "N/A" denotes a convergence failure.
Only the systems that fail to converge or converge to a
total energy above the threshold are given in the tables.

Calculator used: GPAW mode='fd' (see :svn:`gpaw/test/big/scf/g2_1_pbe0_fd.py`)

.. image:: scf_g2_1_pbe0_fd_calculator_steps.png

.. csv-table::
   :file: scf_g2_1_pbe0_fd_energy.csv

DeltaCodesDFT
=============
PBE total energy calculation of bulk systems from
http://molmod.ugent.be/DeltaCodesDFT

All systems are calculated in spin-polarized mode, with the
number of bands specified as ``nbands=-5``, unless stated otherwise.
Corresponding dacapo https://wiki.fysik.dtu.dk/dacapo results given
for comparison.

The runs marked with "initial cg iter N" perform N ``eigensolver=cg``
iterations before switching to the selected eigensolver.
Those initial iterations are included in the final reported number.

Calculator used: GPAW mode=PW() (see :svn:`gpaw/test/big/scf/dcdft_pbe_pw.py`)

.. image:: scf_dcdft_pbe_pw_calculator_steps.png

.. csv-table::
   :file: scf_dcdft_pbe_pw_energy.csv
