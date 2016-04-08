.. module:: ase.calculators.nwchem

======
NWChem
======

`NWChem <http://www.nwchem-sw.org>`_ is a computational chemistry code
based on gaussian basis functions or plane-waves.


Environment variable
====================

.. highlight:: bash

The default command to start NWChem is ``nwchem PREFIX.nw >
PREFIX.out``.  You can change that by setting the environment variable
:envvar:`ASE_NWCHEM_COMMAND`.


Example
=======

Here is an example of how to calculate optimize the geometry of a
water molecule using PBE::

  $ asec H2O optimize -c nwchem -p xc=PBE 
  LBFGS:   0  16:17:29    -2064.914841       1.9673
  LBFGS:   1  16:17:31    -2064.963074       0.9482
  LBFGS:   2  16:17:32    -2064.976603       0.1425
  LBFGS:   3  16:17:33    -2064.977216       0.0823
  LBFGS:   4  16:17:35    -2064.977460       0.0010
  $ ase-gui H2O.traj@-1 -tg "a(1,0,2),d(0,1)"
  102.577881445 1.00806894632
