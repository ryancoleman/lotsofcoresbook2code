Plane wave basis
================

:Who:
    Jens Jørgen

A plane wave implementation is already in trunk.  It's based on FFTW
and does the projector wave function overlaps in reciprocal space with
BLAS's ZGEMM.

Missing features:

* Meta-GGA
* Dipole layer correction
* Parallelization over states works only when number of states is
  divisible by number of processors
* and maybe more ...
