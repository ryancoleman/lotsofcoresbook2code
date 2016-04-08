=======================
xc2.rz.uni-karlsruhe.de
=======================

Here you find information about the the system
`<http://www.rz.uni-karlsruhe.de/ssck/5753.php>`_.

The installation works using ``gcc``::

 module load gcc/4.1.2/default

and numpy-1.0.4 (site.cfg was not used). The installation of gpaw
requires to modify customize.py to::

 libraries = ['acml', 'gfortran']
 library_dirs = ['/software/all/acml/acml4.0/gfortran64/lib']
