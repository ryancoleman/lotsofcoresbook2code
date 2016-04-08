.. _CentOS:

======
CentOS
======

Here you find information about the the system http://www.centos.org/.

.. _PGO_gcc_EL5:

Profile guided optimization
===========================

Example how describes how to use
`profile guided optimization <http://en.wikipedia.org/wiki/Profile-guided_optimization>`_
to compile GPAW with gcc version **4.3.0** on CentOS 5.3:

- starting at :ref:`developer_installation`,
  modify :file:`customize.py` so ``extra_compile_args`` reads::

    opt_string = '-fprofile-generate'
    extra_compile_args =['-std=c99', opt_string]

- moreover, ``mpicompiler`` must be set, and ``mpilinker`` read::

    mpilinker = mpicompiler+' '+opt_string

- build GPAW's :ref:`c_extension` as described at
  :ref:`developer_installation`.  This will create ``*.gcno`` files in
  the :file:`./build/temp.<platform>-<python-version>/c/` directory.

- perform a test run using :file:`gpaw-python`.  This will create
  ``*.gcda`` files in the :file:`./build/temp.<platform>-python-version/c/`
  directory.

- remove object files and :file:`_gpaw.so` (example for **linux-i686** platform, python **2.4**)::

   find build/temp.linux-i686-2.4/ -name "*.o" | xargs rm
   rm -f build/lib.linux-i686-2.4/_gpaw.so

- change :file:`customize.py` so ``opt_string`` reads::

    opt_string = '-fprofile-use'

  and rebuild GPAW's :ref:`c_extension`.


System wide installation with yum
=================================

**Warning**: this section is outdated!

The steps described below require root access and assume bash shell:

- `configure fys yum repository <https://wiki.fysik.dtu.dk/niflheim/Cluster_software_-_RPMS#configure-fys-yum-repository>`_

- on EL/CentOS 6 i386 or x86_64, as root:

  - install gpaw and dependencies::

      yum -y install --enablerepo=fys_el campos-gpaw

  - install optional packages::

      yum -y install --enablerepo=fys_el scipy ScientificPython

- on EL/CentOS 5 i386 or x86_64, as root:

  - install gpaw and dependencies::

      yum -y install --enablerepo=fys_el,epel,atrpms campos-gpaw

  - install optional packages::

      yum -y install --enablerepo=fys_el,epel,atrpms scipy ScientificPython
 
.. note::

   There are no new releases of fys packages after "End of Regular Life Cycle"
   of CentOS releases: see https://access.redhat.com/support/policy/updates/errata/
