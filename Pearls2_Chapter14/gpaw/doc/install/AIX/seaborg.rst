.. _seaborg:

=================
seaborg.nersc.gov
=================

Here you find information about the the system
`<http://www.nersc.gov/nusers/systems/SP/>`_.

We need to use the mpi-enabled compiler ``mpcc`` and we should link to
LAPACK before ESSL.  Make sure LAPACK is added::

  $ module add lapack

and use this customize.py::

  from os import environ
  mpicompiler = 'mpcc'
  libraries = ['f']
  extra_link_args += [environ['LAPACK'], '-lessl']

The Numeric Python extension is not installed on NERSC, so we should
install it.  Get the Numeric-24.2 and do this::

  $ wget http://downloads.sourceforge.net/numpy/Numeric-24.2.tar.gz
  $ gunzip -c Numeric-24.2.tar.gz | tar xf -
  $ cd Numeric-24.2
  $ python setup.py install --home=$HOME

and put the :file:`$HOME/lib/python/Numeric` directory in your
:envvar:`PYTHONPATH`.

Now we are ready to :ref:`compile GPAW <installationguide>`
