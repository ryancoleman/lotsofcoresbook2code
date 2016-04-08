.. _vdwcorrection:

========================
van der Waals correction
========================

A correction on top of the PBE functional has been proposed
by Tkachenko and Scheffler [#TS09]_. While nearly all parameters
are obtained from ab-initio calculations, the method requires
nearly no additional computational cost and performs very well:

============== ===  ===== ====== ======= ====
.              PBE  TPSS  vdW-DF vdW-DF2 TS09
============== ===  ===== ====== ======= ====
Mean deviation 115  154   76     48      15
RMS deviation  108  128   60     42      15
============== ===  ===== ====== ======= ====

Error in energies compared to CCSD results of the S26 test set.
All values in meV.
GPAW calculations were done with h=0.18 and at least 4 A vacuum.
The TS09 results are in good agreement to the results obtained with
the FHI-aims code [#Hanke11jcc]_.

Calculating the S26 test set 
============================

As an example of the usage, here the S26 test set is calculated:

.. literalinclude:: s22_test.py

.. [#TS09] Tkachenko and Scheffler Phys. Rev. Lett. 102 (2009) 073005
.. [#Hanke11jcc] Felix Hanke J. Comp. Chem. 32 (2011) 1424
