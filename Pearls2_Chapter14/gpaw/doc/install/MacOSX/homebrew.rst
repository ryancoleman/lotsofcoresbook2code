.. _homebrew:

========
Homebrew
========

Mountain Lion
=============

Install https://developer.apple.com/xcode/ and activate it from a terminal::

  xcodebuild -license

After installing xcode install also its *Command-Line Tools* (provides
*llvm-gcc compiler* on the command line).
After launching Xcode, in the top menubar, close to the *Apple*, choose
Xcode -> Preferences -> Downloads).

Make sure the compilers are in place::

  which llvm-gcc

Follow the instructions for installing Homebrew http://brew.sh/
the famous::

  ruby -e "$(curl -fsSL https://raw.github.com/Homebrew/homebrew/go/install)"

and configure your init scripts *~/.bash_profile*::

  # Set architecture flags
  export ARCHFLAGS="-arch x86_64"

  # personal installation of pip
  export PATH=/Users/$USER/pip_latest:$PATH
  export PYTHONPATH=/Users/$USER/pip_latest:$PYTHONPATH

  pyver=`python -c "from distutils import sysconfig; print sysconfig.get_python_version()"`

  # pip --user installations of packages
  export PATH=/Users/$USER/Library/Python/${pyver}/bin:$PATH
  export PYTHONPATH=/Users/$USER/Library/Python/${pyver}/lib/python/site-packages:$PYTHONPATH

  # homebrew
  # Ensure user-installed binaries take precedence
  export PATH=/usr/local/share/python:/usr/local/bin:$PATH
  export PYTHONPATH=/usr/local/lib/python${pyver}/site-packages:$PYTHONPATH
  # hack gtk-2.0
  export PYTHONPATH=/usr/local/lib/python${pyver}/site-packages/gtk-2.0:$PYTHONPATH
  # https://github.com/mxcl/homebrew/issues/16891
  export PKG_CONFIG_PATH=`brew --prefix libffi`/lib/pkgconfig:$PKG_CONFIG_PATH
  export PKG_CONFIG_PATH=/opt/X11/lib/pkgconfig:$PKG_CONFIG_PATH
  export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH

  # virtualenv
  # virtualenv should use Distribute instead of legacy setuptools
  export VIRTUALENV_DISTRIBUTE=true
  # Centralized location for new virtual environments
  export PIP_VIRTUALENV_BASE=$HOME/Virtualenvs
  # pip should only run if there is a virtualenv currently activated
  export PIP_REQUIRE_VIRTUALENV=true
  # cache pip-installed packages to avoid re-downloading
  export PIP_DOWNLOAD_CACHE=$HOME/.pip/cache

Verify your homebrew::

  brew doctor

Update with::

  brew update

If you prefer to use OS X python (recommended!), install ``pip``::

  mkdir -p ~/pip_latest
  easy_install --install-dir=$HOME/pip_latest pip  # fails under $HOME/pip

Note that homebrew recommends installing its own python, but by doing so
be prepared on other troubles than addressed in the installation described here.

www.virtualenv.org allows you to run different versions of python modules after
having them configured in different virtualenvs.
It is a convenient way of keeping GPAW with its corresponding
ASE version isolated form the globally installed python modules.

Install virtualenv (you don't need ``--user`` if installing with homebrew python)::

  PIP_REQUIRE_VIRTUALENV=false pip install --user virtualenv
  mkdir ~/Virtualenvs

If you are only installing ASE skip the next section.

Installing GPAW requirements
----------------------------

Install the following homebrew packages::

  brew install gfortran
  brew install gcc
  brew install openmpi
  brew install libxc

Install GPAW setups::

  brew install https://svn.fysik.dtu.dk/projects/gpaw/trunk/doc/install/MacOSX/gpaw-setups.rb

Configure a virtualenv for GPAW, e.g. trunk::

  cd ~/Virtualenvs
  virtualenv gpaw-trunk && cd gpaw-trunk
  . bin/activate

Install the dependencies::

  pip install python-ase
  pip install numpy

and install GPAW (still inside of the virtualenv) accordingly to
:ref:`installationguide_standard` (``python setup.py install``).

Installing ASE requirements
---------------------------

If you prefer to have matplotlib available you need to
install http://xquartz.macosforge.org, reboot, and additionally::

  brew install pygtk

Configure a virtualenv for the latest stable release of ASE::

  cd ~/Virtualenvs
  virtualenv ase && cd ase
  . bin/activate

Now, install ASE inside of virtualenv::

  pip install python-ase
  pip install numpy

Make sure the PKG_CONFIG_PATH correctly
https://github.com/mxcl/homebrew/issues/16891
and then, again inside of virtualenv::

  pip install python-dateutil  # OS X version is outdated!

The latest, precompiled versions of matplotlib (1.3.1) are missing
backend_gdk.so, and therefore compile an older version::

  pip install matplotlib==1.1.1
