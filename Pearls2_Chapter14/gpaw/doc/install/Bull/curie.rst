.. _curie:

=====================================================
curie.ccc.cea.fr (Intel Nehalem, Infiniband QDR, MKL)
=====================================================

Here you find information about the the system
`<http://www-hpc.cea.fr/en/complexe/tgcc-curie.htm>`_.

Numpy is installed system wide (under specific python module), 
so separate installation is not needed.

Building GPAW with gcc
======================

Build GPAW using **gcc** with the configuration file
:svn:`~doc/install/Bull/customize_curie_gcc.py`.

.. literalinclude:: customize_curie_gcc.py

and by executing::

  module unload intel
  module load gnu
  export OMPI_MPICC=gcc
  module load mkl
  module load phdf5
  module load python
  python setup.py install --home=MY_INSTALLATION_DIR

