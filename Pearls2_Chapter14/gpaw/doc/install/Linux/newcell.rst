.. _newcell:

==================
newcell.crc.nd.edu
==================

Here you find information about the the system
`<http://crcmedia.hpcc.nd.edu/wiki/index.php/Available_Hardware>`_.

The installation of user's packages on newcell described below uses
`modules <http://modules.sourceforge.net/>`_, and assumes csh shell:

- packages are installed under ``~/CAMd``::

   mkdir ~/CAMd
   cd ~/CAMd

- module files are located under ``~/CAMd/modulefiles``

- download the :svn:`~doc/install/Linux/customize_newcell.py` file::

   wget https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Linux/customize_newcell.py

  .. literalinclude:: customize_newcell.py

- install packages, deploy modules and test with :svn:`~doc/install/Linux/install_newcell.csh`,
  buy running ``csh install_newcell.csh``:

  .. literalinclude:: install_newcell.csh

  **Note** that every time you wish to install a new version of a package,
  and deploy new module file, better keep the old module file.

- submit the test job::

   qsub submit.csh

  using the following :file:`submit.csh`::

   #!/bin/csh

   #$ -q short
   #$ -pe ompi-8 8
   #$ -l arch=lx24-amd64
   ##$ -M jbray2@nd.edu
   ##$ -m abe

   module load ompi/1.3.2-gnu

   module use --append /afs/crc.nd.edu/user/j/jbray2/CAMd/modulefiles
   module load nose
   module load numpy
   module load campos-ase3
   module load campos-gpaw-setups
   module load campos-gpaw

   setenv P4_GLOBMEMSIZE 268435456
   setenv P4_SOCKBUFSIZE 262144

   echo job was accepted on:
   date

   mpirun -np $NSLOTS `which gpaw-python` `which gpaw-test`

   echo Job has completed.
   date

   exit 0
