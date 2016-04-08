.. _exx:

==============
Exact exchange
==============

**THIS PAGE IS OUTDATED**

Inclusion of the non-local Fock operator as an exchange-correclation functional is an experimental feature in gpaw.

The current implementation *lacks* the following features:

* Support for periodic systems.
   Actually, the code won't complain, but the results have not been tested.
* Support for k-point sampling.
   No consideration has been made as to multiple k-points, or even comlplex wave functions, so this definitely won't work.
* Forces.
   Force evaluations when including (a fraction of -) the fock operator in the xc-functional has been implemented, but has not been tested.
* Fractional occupations.
   For technical reasons, we had to decouple occupied and unoccupied states. This makes fractional occupations imposible. (Warning: the code will not raise an exception, but probably won't converge).
* Speed.
   Inclusion of Fock exchange is exceedingly slow. The bottleneck is solving the poisson integrals of the Fock operator, which is currently done using an iterative real-space solver with a zero initial guess for the potential at each SCF cycle. This chould be optimized.

One way to speed up an exact-exchange (or hybrid) calculation is to use the coarse grid (used for wave functions) instead of the finegrid (used for for densities) for the Fock potentials. This should give a speed-up factor of 8. This can be specified in the ``xc`` keyword like in this example :svn:`~gpaw/test/exx_coarse.py`

Parallelization using domain decomposition is fully supported.

The Fock operator can be used to do the hybrid functional PBE0, and of course, Hartree-Fock type EXX. These are accessed by setting the ``xc`` keyword to ``PBE0`` or ``EXX`` respectively.

A thesis on the implementation of EXX in the PAW framework, and the
specifics of the GPAW project can be seen on the :ref:`literature
<literature_reports_presentations_and_theses>` page.

A comparison of the atomization energies of the g2-1 test-set calculated in VASP, Gaussian03, and GPAW is shown in the below two figures for the PBE and the PBE0 functional respectively.

.. image:: g2test_pbe.png

.. image:: g2test_pbe0.png

In the last figure, the curve marked ``GPAW (nonself.)`` is a non-selfconsistent PBE0 calculation using self-consistent PBE orbitals.
