.. _odyssey:

=======================
odyssey.fas.harvard.edu
=======================

Information about the Harvard Odyssey cluster can be found at
`<http://rc.fas.harvard.edu>`_.

Intel C++ compiler, MKL, and OpenMPI are required::

> module load hpc/numpy-1.6.0_python-2.7.1
> module load hpc/intel-and-mkl-12.3.174
> module load hpc/openmpi-1.5.3_intel-12.3.174ib

The following :file:`customize.py` file can be used to build GPAW:

.. literalinclude:: customize_odyssey.py



