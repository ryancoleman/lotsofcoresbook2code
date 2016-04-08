.. _dcdft:

===============
Delta Codes DFT
===============

Codes precision estimated for PBE exchange-correlation functional
on the database of bulk systems from http://molmod.ugent.be/DeltaCodesDFT.

The Delta values calculated based on the EOS parameters:
*V0* in Å**3/atom, bulk modulus (*B0*) in GPa, and
pressure derivative of the bulk modulus *B1* (dimensionless).
Percentage errors with respect to http://www.wien2k.at/.

**Note**: in order to get an impression about the magnitude of the
Delta factors one could for example divide them by the corresponding
Deltas of PBE exchange-correlation functional with respect to
experiment, available in the publication
http://molmod.ugent.be/DeltaCodesDFT.  This can be
performed using the ``calcDelta.py`` script provided at the website::

  python calcDelta.py VASP-relaxed.txt exp.txt --stdout | grep -E -v "#|-|np" > test.txt
  python calcDelta.py VASP.txt WIEN2k.txt --stdout | grep -E -v "#|-|np" | join test.txt - | grep -v "N/A" | awk '{print $1, $3/$2*100}'


GPAW PAW
--------

Calculated with: :svn:`gpaw/test/big/dcdft/pbe_gpaw_pw.py`.

EOS
+++

.. csv-table::
   :file: dcdft_pbe_gpaw_pw_raw.csv

Delta precision measure
+++++++++++++++++++++++

Consult http://molmod.ugent.be/DeltaCodesDFT


FHI-AIMS tight
--------------

Calculated with: :svn:`gpaw/test/big/dcdft/pbe_aims.py`.

EOS
+++

.. csv-table::
   :file: dcdft_aims.tight.01.16.db_raw.csv

Delta precision measure
+++++++++++++++++++++++

Consult http://molmod.ugent.be/DeltaCodesDFT


Dacapo
------

https://wiki.fysik.dtu.dk/dacapo
Calculated with: :svn:`gpaw/test/big/dcdft/pbe_jacapo.py`.

EOS
+++

.. csv-table::
   :file: dcdft_pbe_jacapo_raw.csv

Delta precision measure
+++++++++++++++++++++++

Consult http://molmod.ugent.be/DeltaCodesDFT


Abinit FHI
----------

Abinit 5.4.4p, `GGA_FHI <http://www.abinit.org/downloads/psp-links/gga_fhi>`_
pseudopotentials, calculated with: :svn:`gpaw/test/big/dcdft/pbe_abinit_fhi.py`.

EOS
+++

.. csv-table::
   :file: dcdft_pbe_abinit_fhi_raw.csv

Delta precision measure
+++++++++++++++++++++++

Consult http://molmod.ugent.be/DeltaCodesDFT
