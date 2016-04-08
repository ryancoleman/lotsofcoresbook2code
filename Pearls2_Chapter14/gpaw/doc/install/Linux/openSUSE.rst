.. _openSUSE:

========
openSUSE
========

Here you find information about the the system
`<http://www.opensuse.org/>`_.

System wide installation with yast
==================================

**Warning**: this section is outdated!

The steps described below require root access and assume bash shell:

- `configure fys yum repository <https://wiki.fysik.dtu.dk/niflheim/Cluster_software_-_RPMS#configure-fys-yum-repository>`_

- on openSUSE 11.2 or newer i386 or x86_64, as root:

  - list available packages::

      zypper packages -r fys_openSUSE

  - install gpaw and dependencies::

      yast -i campos-gpaw

  - install optional packages::

      yast -i scipy ScientificPython

.. note::

   There are no new releases of fys packages for "Discontinued distributions"
   of openSUSE: see http://en.opensuse.org/Lifetime

