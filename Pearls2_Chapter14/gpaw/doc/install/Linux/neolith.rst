.. _neolith:

==================
neolith.nsc.liu.se
==================

Here you find information about the the system
http://www.nsc.liu.se/systems/neolith.

The installation of user's packages on neolith described below uses
`cmod <http://www.lysator.liu.se/cmod/>`_ modules:

- packages are installed under ``~/apps``::

   mkdir ~/apps

- module files are located under ``~/modulefiles``::

   mkdir ~/modulefiles

- build `nose <http://code.google.com/p/python-nose/>`_::

   cd ~/apps
   wget http://python-nose.googlecode.com/files/nose-0.10.1.tar.gz
   tar zxf nose-0.10.1.tar.gz
   cd nose-0.10.1
   python setup.py install --root=~/apps/nose-0.10.1-1

- build `numpy <http://numpy.scipy.org/>`_::

   cd ~/apps
   wget http://downloads.sourceforge.net/numpy/numpy-1.3.0.tar.gz
   tar zxf numpy-1.3.0.tar.gz
   cd  numpy-1.3.0
   python setup.py install --root=~/apps/numpy-1.3.0-1

- build `ase <https://wiki.fysik.dtu.dk/ase/>`_::

   cd ~/apps
   svn co https://svn.fysik.dtu.dk/projects/ase/trunk ase
   cd ase
   python setup.py sdist; cp dist/python-ase-*.tar.gz ..
   cd ..
   tar zxf python-ase-3.2.0.962.tar.gz

- build gpaw-setups::

   cd ~/apps
   wget --no-check-certificate "http://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.5.3574.tar.gz"
   tar zxf gpaw-setups-0.5.3574.tar.gz

- deploy modules::

    export MODULEFILES="/home/x_andke/modulefiles"
    export APPS="/home/x_andke/apps"

    mkdir ${MODULEFILES}/numpy
    cat <<EOF > ${MODULEFILES}/numpy/1.3.0-1.el5.gfortran.python2.4.default.blas.lapack
    \$(/etc/cmod/modulegroups numpy numpy/1.3.0-1.el5.gfortran.python2.4.default.blas.lapack)
    #prereq python-nose # nose
    prepend-path    PATH            ${APPS}/numpy-1.3.0-1/usr/bin
    prepend-path    PYTHONPATH      ${APPS}/numpy-1.3.0-1/usr/lib64/python2.4/site-packages
    EOF
    ln -s ${MODULEFILES}/numpy/1.3.0-1.el5.gfortran.python2.4.default.blas.lapack ${MODULEFILES}/numpy/default

    mkdir ${MODULEFILES}/nose
    cat <<EOF > ${MODULEFILES}/nose/0.10.1-1.el5.gfortran.python2.4
    \$(/etc/cmod/modulegroups nose nose/0.10.1-1.el5.gfortran.python2.4)
    prepend-path    PATH            ${APPS}/nose-0.10.1-1/usr/bin
    prepend-path    PYTHONPATH      ${APPS}/nose-0.10.1-1/usr/lib/python2.4/site-packages
    EOF
    ln -s ${MODULEFILES}/nose/0.10.1-1.el5.gfortran.python2.4 ${MODULEFILES}/nose/default

    mkdir ${MODULEFILES}/campos-ase3
    cat <<EOF > ${MODULEFILES}/campos-ase3/3.2.0.962-1.el5.python2.4
    \$(/etc/cmod/modulegroups campos-ase3 campos-ase3/3.2.0.962-1.el5.python2.4)
    #prereq numpy # numpy
    prepend-path    PATH            ${APPS}/python-ase-3.2.0.962/tools
    prepend-path    PYTHONPATH      ${APPS}/python-ase-3.2.0.962
    EOF
    ln -s ${MODULEFILES}/campos-ase3/3.2.0.962-1.el5.python2.4 ${MODULEFILES}/campos-ase3/default

    mkdir ${MODULEFILES}/campos-gpaw-setups
    cat <<EOF > ${MODULEFILES}/campos-gpaw-setups/0.5.3574-1.el5
    \$(/etc/cmod/modulegroups campos-gpaw-setups campos-gpaw-setups/0.5.3574-1.el5)
    prepend-path    GPAW_SETUP_PATH      ${APPS}/gpaw-setups-0.5.3574
    EOF
    ln -s ${MODULEFILES}/campos-gpaw-setups/0.5.3574-1.el5 ${MODULEFILES}/campos-gpaw-setups/default

    mkdir ${MODULEFILES}/campos-gpaw
    cat <<EOF > ${MODULEFILES}/campos-gpaw/0.6.3934-1.el5.gfortran.python2.4.openmpi.mkl.10.0.4.023.mkl_lapack
    \$(/etc/cmod/modulegroups campos-gpaw campos-gpaw/0.6.3934-1.el5.gfortran.python2.5.openmpi.mkl.10.0.4.023.mkl_lapack)
    #prereq numpy # numpy
    #prereq campos-ase3 # ase
    #prereq campos-gpaw-setups # gpaw-setups
    prepend-path    PATH            ${APPS}/gpaw-0.6.3934/tools
    prepend-path    PATH            ${APPS}/gpaw-0.6.3934/build/bin.linux-x86_64-2.4
    prepend-path    PYTHONPATH      ${APPS}/gpaw-0.6.3934
    EOF
    ln -s ${MODULEFILES}/campos-gpaw/0.6.3934-1.el5.gfortran.python2.4.openmpi.mkl.10.0.4.023.mkl_lapack ${MODULEFILES}/campos-gpaw/default

  **Note** that every time you wish to install a new version of a package,
  you (``svn up`` is necessary), create new tarball,
  and deploy new module file, keeping the old module file.

- test numpy installation::

   module use ${MODULEFILES}
   module load nose
   module load numpy
   python -c "import numpy; numpy.test()"

- use :svn:`~doc/install/Linux/customize_neolith.py`:

  .. literalinclude:: customize_neolith.py

  to build `gpaw <https://wiki.fysik.dtu.dk/gpaw/>`_::

   cd ~/apps
   svn co https://svn.fysik.dtu.dk/projects/gpaw/trunk gpaw
   cd gpaw
   python setup.py sdist; cp dist/gpaw-*.tar.gz ..
   cd ..
   tar zxf gpaw-0.6.3934.tar.gz
   cd gpaw-0.6.3934
   wget https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Linux/customize_neolith.py -O customize.py
   module load openmpi/1.2.7-i101017
   python setup.py build_ext --remove-default-flags

- test gpaw installation by loading the modules::

   module load campos-ase3
   module load campos-gpaw-setups
   module load campos-gpaw
   export OMP_NUM_THREADS=1

  and :ref:`running_tests`.

- **logout**, and login again.

- submit a test job::

   cp ~/apps/gpaw-0.6.3934/test/CH4.py ~/
   cd
   sbatch -N 1 --tasks-per-node 4 submit.sh

  using the following :file:`submit.sh`::

   #!/bin/bash
   #SBATCH -N 1
   #SBATCH -t 00:10:00

   export OMP_NUM_THREADS=1
   . /etc/cmod/path.sh
   module use /home/x_andke/modulefiles
   module load openmpi/1.2.7-i101017
   module load nose
   module load numpy
   module load campos-ase3
   module load campos-gpaw-setups
   module load campos-gpaw

   mpprun --force-mpi="openmpi/1.2.7-i101017" `which gpaw-python` ./CH4.py --sl_diagonalize=2,1,2

