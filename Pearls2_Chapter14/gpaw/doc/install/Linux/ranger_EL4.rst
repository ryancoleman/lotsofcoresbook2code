.. _ranger:

======================
ranger.tacc.utexas.edu
======================

Here you find information about the the system
http://services.tacc.utexas.edu/index.php/ranger-user-guide.

The installation of user's packages on ranger EL4, 64-bit described below uses
`modules <http://modules.sourceforge.net/>`_, and assumes *csh* shell:

- packages are installed under ``~/CAMd``::

   mkdir ~/CAMd
   cd ~/CAMd

- module files are located under ``~/CAMd/modulefiles``

- download the :svn:`~doc/install/Linux/customize_ranger_EL4.py` file::

   wget https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Linux/customize_ranger_EL4.py

  .. literalinclude:: customize_ranger_EL4.py

- download packages with :svn:`~doc/install/Linux/download_ranger.sh`,
  buy running ``sh download_ranger.sh``:

  .. literalinclude:: download_ranger.sh

- start with the following modules loaded::

   login3% module list
   Currently Loaded Modules:
     1) TACC-paths           9) srb-client/3.4.1  17) globus/4.0.8
     2) Linux                10) tg-policy/0.2    18) GLOBUS-4.0
     3) cluster-paths        11) tgproxy/0.9.1    19) TERAGRID-DEV
     4) pgi/7.2-5            12) tgresid/2.3.4    20) CTSSV4
     5) mvapich/1.0.1        13) tgusage/3.0      21) gzip/1.3.12
     6) binutils-amd/070220  14) uberftp/2.4      22) tar/1.22
     7) TERAGRID-paths       15) tginfo/1.0.1     23) cluster
     8) gx-map/0.5.3.3       16) TERAGRID-BASIC   24) TACC

- unload/load the modules::

    module switch pgi/7.2-5 gcc/4.4.5
    module swap mvapich openmpi/1.3b
    module load python/2.5.2
    module load mkl/10.0

- install packages, deploy modules and test with :svn:`~doc/install/Linux/install_ranger_EL4.sh`,
  buy running ``sh install_ranger_EL4.sh``:

  .. literalinclude:: install_ranger_EL4.sh

  **Note** that every time you wish to install a new version of a package,
  and deploy new module file, better keep the old module file.

- submit the test job::

   qsub submit.sh

  using the following :file:`submit.sh`::

   #!/bin/bash      

   #$ -V   # Inherit the submission environment
   #$ -cwd         # Start job in submission directory
   ##$ -N myMPI    # Job Name
   ##$ -j y        # Combine stderr and stdout
   ##$ -o $JOB_NAME.o$JOB_ID       # Name of the output file (eg. myMPI.oJobID)
   #$ -pe 16way 32         # Requests 16 tasks/node, 32 cores total
   #$ -q development       # OR Queue name "normal"
   #$ -l h_rt=00:40:00     # Run time (hh:mm:ss) - 40 mins
   ##$ -M  # Use email notification address
   ##$ -m be       # Email at Begin and End of job
   #set -x         # Echo commands, use "set echo" with csh

   module use --append /share/home/01067/tg803307/CAMd/modulefiles
   module load python/2.5.2
   module load nose/0.11.3-1
   module load numpy/1.5.0-1
   module load campos-ase3
   module load campos-gpaw-setups
   module unload pgi
   module load gcc/4.4.5
   module unload mvapich
   module load openmpi/1.3b
   module load mkl/10.0
   module load campos-gpaw

   # wget http://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/devel/256H2O/b256H2O.py

   ibrun `which gpaw-python` b256H2O.py

- to enable the installation permanently add the following to *~/.bashrc*::

   module use --append /share/home/01067/tg803307/CAMd/modulefiles
   module load python/2.5.2
   module load nose/0.11.3-1
   module load numpy/1.5.0-1
   module load campos-ase3
   module load campos-gpaw-setups
   module unload pgi
   module load gcc/4.4.5
   module unload mvapich
   module load openmpi/1.3b
   module load mkl/10.0
   module load campos-gpaw
