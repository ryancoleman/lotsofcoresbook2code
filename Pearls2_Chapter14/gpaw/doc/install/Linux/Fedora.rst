.. _Fedora:

======
Fedora
======

Here you find information about the the system
`<http://fedoraproject.org/>`_.

System wide installation with yum
=================================

**Warning**: this section is outdated!

The steps described below require root access and assume bash shell:

- `configure fys yum repository <https://wiki.fysik.dtu.dk/niflheim/Cluster_software_-_RPMS#configure-fys-yum-repository>`_

- on Fedora 12 or newer i386 or x86_64, as root:

  - install gpaw and dependencies::

      yum -y install --enablerepo=fys_fc campos-gpaw

  - install optional packages::

      yum -y install scipy ScientificPython

.. note::

   There are no new releases of fys packages for "Old Unsupported Releases"
   of Fedora: see http://fedoraproject.org/wiki/Releases
 
