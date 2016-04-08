.. _jump:

==================
jump.fz-juelich.de
==================

Here you find information about the the system
`<http://www.fz-juelich.de/jsc/service/sco_ibmP6>`_.

The only way we found to compile numpy is using python2.3 and
numpy-1.0.4. The next version numpy-1.1.0 did not work
unfortunately. In addition the usage of the generic IBM lapack/blas in
numpy does not work, hence one has to use site.cfg::

  : diff site.cfg site.cfg.example
  58,60c58,60
  < [DEFAULT]
  < library_dirs =
  < include_dirs =
  ---
  > #[DEFAULT]
  > #library_dirs = /usr/local/lib
  > #include_dirs = /usr/local/include

With his change numpy compiles and the installation of ASE and gpaw
does not cause problems.

