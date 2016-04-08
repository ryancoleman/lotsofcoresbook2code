.. _FreeBSD:

=======
FreeBSD
=======

Here you find information about the the system
`<http://freebsd.org/>`_.

To build gpaw add to the ``gpaw/customize.py``::

  compiler='gcc44'
  extra_compile_args += ['-Wall', '-std=c99'] 
  library_dirs += ['/usr/local/lib']
  libraries += ['blas', 'lapack', 'gfortran']

If you want to build a parallel version install mpich2 (net/mpich2). Openmpi
does currently not work with gpaw!
