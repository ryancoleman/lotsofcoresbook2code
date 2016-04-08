.. _Ubuntupackage:

==============
Ubuntu package
==============

Here you find information about the system: `<http://www.ubuntu.com/>`_.

GPAW is available as an `Ubuntu package
<https://launchpad.net/~campos-dev/+archive/campos>`_ on `Launchpad
<https://launchpad.net/>`_ for Ubuntu 9.10 or newer. To install:

- Add the package archive to the system's software
  sources::

    sudo add-apt-repository ppa:campos-dev/campos

  .. note::

    More up-to-date packages can be usually found at::

      sudo add-apt-repository ppa:askhl/ppa

- Update the package cache::

    sudo apt-get update

- Install GPAW::

    sudo apt-get install gpaw

This will also install ASE, the GPAW setups, MPI and other
dependencies.  ASE and GPAW can be automatically upgraded when a new
stable version is released.


