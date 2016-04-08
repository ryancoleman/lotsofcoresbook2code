.. _nanolab:

=======================
nanolab.cnf.cornell.edu
=======================

Here you find information about the the system
http://www.cnf.cornell.edu/cnf5_tool.taf?_function=detail&eq_id=111.

The installation of user's packages on nanolab EL4, 32-bit described below uses
`modules <http://modules.sourceforge.net/>`_, and assumes *bash* shell:

- packages are installed under ``~/CAMd``::

   mkdir ~/CAMd
   cd ~/CAMd

- module files are located under ``~/CAMd/modulefiles``

- download the :svn:`~doc/install/Linux/customize_nanolab_EL4.py` file::

   wget https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Linux/customize_nanolab_EL4.py

  .. literalinclude:: customize_nanolab_EL4.py

- download packages with :svn:`~doc/install/Linux/download_nanolab.sh`,
  buy running ``sh download_nanolab.sh``:

  .. literalinclude:: download_nanolab.sh

- from *nanolab.cnf.cornell.edu* login to one of c-nodes (Red Hat 4, 32-bit)::

    ssh c7.cnf.cornell.edu

- install packages, deploy modules and test with
  :svn:`~doc/install/Linux/install_nanolab_EL4.sh`, buy running ``sh
  install_nanolab_EL4.sh``:

  .. literalinclude:: install_nanolab_EL4.sh

  **Note** that every time you wish to install a new version of a package,
  and deploy new module file, better keep the old module file.


- submit the test job::

   qsub submit.sh

  using the following :file:`submit.sh`::

   TODO

- to enable the installation permanently add the following to *~/.bashrc*::

   module use --append /home/karsten/CAMd/modulefiles
   module load numpy
   module load campos-ase3
   module load campos-gpaw-setups
   module load intel_compilers/11.1
   module load openmpi/1.3.3
   module load campos-gpaw
