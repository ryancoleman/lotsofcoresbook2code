.. _juropatest:

======================
juropatest @ FZ-Jülich
======================

Here you find information about the the system
http://www.fz-juelich.de/ias/jsc/juropatest

Pre-installed versions
======================

You may use the pre-installed versions::

  module load intel-para
  module load GPAW

In case you are happy with these versions, you need to install
the setups (next point) and you are done.

Setups
======

The setups are not defined in the per-installed vesrion, so we need
to install them ourselves::

  cd
  GPAW_SETUP_SOURCE=$PWD/source/gpaw-setups
  mkdir -p $GPAW_SETUP_SOURCE
  cd $GPAW_SETUP_SOURCE
  wget https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.9.11271.tar.gz
  tar xzf gpaw-setups-0.9.11271.tar.gz
  
Let gpaw know about the setups::
  
  export GPAW_SETUP_PATH=$GPAW_SETUP_SOURCE/gpaw-setups-0.9.11271

Using the module environment
============================

It is very handy to add our installation to the module environment::

  cd
  mkdir -p modulefiles/gpaw-setups
  cd modulefiles/gpaw-setups
  echo -e "#%Module1.0\nprepend-path       GPAW_SETUP_PATH    $GPAW_SETUP_SOURCE/gpaw-setups-0.9.11271" > 0.9.11271
  
We need to let the system know about our modules::

  module use $HOME/modulefiles

such that we also see them with::

  module avail

Building from trunk
===================

In case that you need a newer version than is installed you might want 
to install gpaw yourself.

We first create a place for gpaw and get the trunk version::

  cd
  GPAW_SOURCE=$PWD/source/gpaw
  mkdir -p $GPAW_SOURCE
  cd $GPAW_SOURCE
  svn checkout https://svn.fysik.dtu.dk/projects/gpaw/trunk trunk

The current trunk version can then be updated by::

  cd $GPAW_SOURCE/trunk
  svn up

We use the installed versions of ASE and libxc::

  module load intel-para
  module load ASE
  module load libxc

and install using
:svn:`~doc/install/Linux/customize_juropatest.py`::

  cd $GPAW_SOURCE/trunk
  mkdir install
  cp customize_juropatest.py customize.py
  python setup.py install --prefix=$PWD/install

It is very handy to add our installation to the module environment::

  cd
  mkdir -p modulefiles/gpaw
  
and the module file  :file:`trunk` should read::

  #%Module1.0

  if {![is-loaded intel-para]} {module load intel-para}
  if {![is-loaded ASE]} {module load ASE}
  if {![is-loaded libxc]} {module load libxc}
  if {![is-loaded gpaw-setups]}  {module load gpaw-setups}

  # change this to your $HOME
  set HOME /homea/hfr08/hfr080

  set gpawhome $HOME/source/gpaw/trunk/install
  prepend-path    PATH                 $gpawhome/bin
  prepend-path    PYTHONPATH           $gpawhome/lib/python2.7/site-packages/
  setenv          GPAW_PYTHON          $gpawhome/bin/gpaw-python

Execution
=========

Job scripts can be written using::

  gpaw-runscript -h

