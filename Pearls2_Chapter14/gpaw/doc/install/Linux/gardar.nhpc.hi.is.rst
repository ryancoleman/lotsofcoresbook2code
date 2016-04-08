.. _gardar.nhpc.hi.is:

=================
gardar.nhpc.hi.is
=================

Information about the machine http://nhpc.hi.is/content/nhpc-system

Instructions assume **bash**, installation under *${HOME}/global/apps*.

Make sure the following modules are loaded::

  $ module list

  Currently Loaded Modules:
    1) intel/13.1   3) openmpi/1.6.5   5) gold/2.2.0.5
    2) mkl/11.0.0   4) python/2.7.5

Setup the root directory::

  export APPHOME=${HOME}/global/apps

  mkdir -p ${APPHOME}
  cd ${APPHOME}

  export GPAW_PLATFORM=`python -c "from distutils import util, sysconfig; print util.get_platform()+'-'+sysconfig.get_python_version()"`

Download software::

  svn co -r 3906 https://svn.fysik.dtu.dk/projects/ase/trunk ase.3906
  svn co -r 12224 https://svn.fysik.dtu.dk/projects/gpaw/trunk gpaw.12224
  wget https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.9.11271.tar.gz
  wget https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.9.9672.tar.gz
  tar zxf gpaw-setups-0.9.11271.tar.gz
  tar zxf gpaw-setups-0.9.9672.tar.gz
  wget "http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.2.1.tar.gz" -O libxc-2.2.1.tar.gz
  tar zxf libxc-2.2.1.tar.gz

and create the modules::

  mkdir -p $APPHOME/modulefiles/{ase,gpaw,gpaw-setups,libxc}

  cat <<EOF > modulefiles/gpaw/0.11.0.12224
  set modulefile [lrange [split [module-info name] {/}] 0 0]
  set release    [lrange [split [module-info name] {/}] 1 1]
  set apphome    $APPHOME/gpaw.12224
  set appname    "gpaw"

  prepend-path PATH \$apphome/build/bin.$GPAW_PLATFORM
  prepend-path PATH \$apphome/tools
  prepend-path PYTHONPATH \$apphome
  prepend-path PYTHONPATH \$apphome/build/lib.$GPAW_PLATFORM
  EOF

  cat <<EOF > modulefiles/ase/3.9.0.3906
  set modulefile [lrange [split [module-info name] {/}] 0 0]
  set release    [lrange [split [module-info name] {/}] 1 1]
  set apphome    $APPHOME/ase.3906
  set appname    "ase"

  prepend-path PATH \$apphome/tools
  prepend-path PYTHONPATH \$apphome
  EOF

  cat <<EOF > modulefiles/gpaw-setups/0.9.9672
  set modulefile [lrange [split [module-info name] {/}] 0 0]
  set release    [lrange [split [module-info name] {/}] 1 1]
  set apphome    $APPHOME/gpaw-setups-0.9.9672
  set appname    "gpaw-setups"

  prepend-path GPAW_SETUP_PATH \$apphome
  EOF

  cat <<EOF > modulefiles/gpaw-setups/0.9.11271
  set modulefile [lrange [split [module-info name] {/}] 0 0]
  set release    [lrange [split [module-info name] {/}] 1 1]
  set apphome    $APPHOME/gpaw-setups-0.9.11271
  set appname    "gpaw-setups"

  prepend-path GPAW_SETUP_PATH \$apphome
  EOF

  cat <<EOF > modulefiles/libxc/2.2.1-1
  set modulefile [lrange [split [module-info name] {/}] 0 0]
  set release    [lrange [split [module-info name] {/}] 1 1]
  set apphome    $APPHOME/libxc-2.2.1-1
  set appname    "libxc"

  prepend-path    LD_LIBRARY_PATH    \$apphome/lib
  prepend-path    PATH               \$apphome/bin
  prepend-path    C_INCLUDE_PATH     \$apphome/include
  prepend-path    PKG_CONFIG_PATH    \$apphome/lib/pkgconfig
  EOF

Build libxc::

  cd $APPHOME/libxc-2.2.1
  ./configure --prefix $APPHOME/libxc-2.2.1-1 --enable-shared
  make
  make install

and GPAW::

  cd $APPHOME
  wget https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/Linux/customize_gardar.py
  cd $APPHOME/gpaw.12224
  python setup.py --remove-default-flags --customize=../customize_gardar.py build_ext 2>&1 | tee build_ext.log

The :file:`customize_gardar.py` looks like:

.. literalinclude:: customize_gardar.py

GPAW tests :file:`gpaw-test` can be submitted like this::

  qsub -l nodes=1:ppn=8 run.sh

where :file:`run.sh` looks like this::

  #!/bin/sh

  . /opt/lmod/lmod/init/sh
  module use --append ~/global/apps/modulefiles
  module load ase/3.9.0.3906
  module load gpaw-setups/0.9.11271
  module load libxc/2.2.1-1
  module load gpaw/0.11.0.12224
  export OMP_NUM_THREADS=1

  mpiexec `which gpaw-python` `which gpaw-test`
