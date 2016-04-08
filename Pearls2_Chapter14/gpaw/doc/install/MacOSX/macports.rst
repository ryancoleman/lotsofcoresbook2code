.. _macports:

========
MacPorts
========

Snow Leopard
============

Follow the instructions here for installing MacPorts and all
supporting packages:
https://sites.google.com/site/naromero/code-development-on-mac-os-x.

This build of GPAW will use MPICH2 1.2, Python 2.7.x and NumPy 1.5.x.

Use the following customize file :svn:`customize_snowleopard_macports.py` to

.. literalinclude:: customize_snowleopard_macports.py

There have been some problems with vecLib library which is linked by
default. Disabling this is in ``config.py`` is straigtforward.  Here
is a diff :svn:`config.disable_vecLib.py.diff`

.. literalinclude:: config.disable_vecLib.py.diff

Then build::

  python setup.py build_ext --customize=customize_snowleopard_macports.py

For Apple users, the MacPorts Project provides a straight-forward
route to obtain all necessary requirements. Unfortunately, MacPorts
does not install the *gtk* bindings to matplotlib by default, which
are required to open the GUI. To get all the ASE prerequisites for
python 2.7 in one single command anyway, install MacPorts and then run::

  $ port install py27-matplotlib +gtk2
