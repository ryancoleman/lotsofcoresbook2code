.. _batman:

==================
batman.chem.jyu.fi
==================

To prepare the compilation, we need to load the required modules and
clean the environment::

 > module purge # remove all modules
 > module add mpt
 > module add mkl
 > unset CC CFLAGS LDFLAGS

We have to change :file:`customize.py` to get the libs and the right compiler::

 # uncomment and change in customize.py
 libraries += ['mpi','mkl']

 mpicompiler = 'gcc'
 custom_interpreter = True

Then compile as usual (``python setup.py build``). This will build the
custom python interpreter for parallel use also.

