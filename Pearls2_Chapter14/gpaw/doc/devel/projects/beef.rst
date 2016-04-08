Bayesian error estimation functional (BEEF)
===========================================

:Who:
    Jess, Keld

We use methods inspired from Bayesian statistics and a large number of datasets
of molecular, surface chemical, and solid state materials properties
to construct new exchange-correlation density functionals.
An essential feature of these functionals is an ensemble of functionals
around the optimum one, which allows an estimate of the computational error
to be easily calculated in a non-self-consistent fashion.

We have recently designed the BEEF-vdW, a GGA with vdW-DF2 type nonlocal
correlation. It is fully implemented in GPAW,
as is the ensemble error estimate.

Presently and in the future we will expand on this work by considering:

* metaGGA density functionals, which us the klinetic energy density
* self-interaction correction methods, e.g., +U, non-Koopman corrections, etc.

Most of the GPAW code related to this project is in
:svn:`~gpaw/gpaw/xc/bee.py` and :svn:`~gpaw/c/ensemble_gga.c`.
