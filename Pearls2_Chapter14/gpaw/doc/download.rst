.. _download:

========
Download
========

.. note::

   GPAW requires ASE. :ase:`Install ASE <download.html>`.

   When done, to determine which way of installing GPAW suits you best
   please read carefully :ref:`installationguide` first!

.. _latest_stable_release:

Latest stable release
=====================

The latest stable release can be obtained from ``svn`` or as a ``tarball``.

===========  =======  ========  =========================  ====================
Date         Version  Revision  Tarfile                    Required ASE version
===========  =======  ========  =========================  ====================
Apr  8 2014  0.10.0_  11364     gpaw-0.10.0.11364.tar.gz_  3.8.1
Mar  7 2012  0.9.0_   8965      gpaw-0.9.0.8965.tar.gz_    3.6.0
May 25 2011  0.8.0_   8092      gpaw-0.8.0.8092.tar.gz_    3.5.1
Aug 11 2010  0.7.2_   6974      gpaw-0.7.2.6974.tar.gz_    3.4.1
Apr 23 2010  0.7_     6383      gpaw-0.7.6383.tar.gz_      3.4.0
Oct  9 2009  0.6_     5147      gpaw-0.6.5147.tar.gz_      3.2.0
Apr  1 2009  0.5_     3667      gpaw-0.5.3667.tar.gz_      3.1.0
Nov 16 2008  0.4_     2734      gpaw-0.4.2734.tar.gz_      3.0.0
===========  =======  ========  =========================  ====================

To check out the latest stable version from SVN, do this::

  $ svn co -r 11364 https://svn.fysik.dtu.dk/projects/gpaw/tags/0.10.0 gpaw-0.10.0

.. _0.10.0:
    https://trac.fysik.dtu.dk/projects/gpaw/browser/tags/0.10.0

.. _gpaw-0.10.0.11364.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-0.10.0.11364.tar.gz

.. _0.9.0:
    https://trac.fysik.dtu.dk/projects/gpaw/browser/tags/0.9.0

.. _gpaw-0.9.0.8965.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-0.9.0.8965.tar.gz

.. _0.8.0:
    https://trac.fysik.dtu.dk/projects/gpaw/browser/tags/0.8.0

.. _gpaw-0.8.0.8092.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-0.8.0.8092.tar.gz

.. _0.7.2:
    https://trac.fysik.dtu.dk/projects/gpaw/browser/tags/0.7.2

.. _gpaw-0.7.2.6974.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-0.7.2.6974.tar.gz

.. _0.7:
    https://trac.fysik.dtu.dk/projects/gpaw/browser/tags/0.7

.. _gpaw-0.7.6383.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-0.7.6383.tar.gz

.. _0.6:
    https://trac.fysik.dtu.dk/projects/gpaw/browser/tags/0.6

.. _gpaw-0.6.5147.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-0.6.5147.tar.gz

.. _0.5:
    https://trac.fysik.dtu.dk/projects/gpaw/browser/tags/0.5

.. _gpaw-0.5.3667.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-0.5.3667.tar.gz

.. _0.4:
    https://trac.fysik.dtu.dk/projects/gpaw/browser/tags/0.4

.. _gpaw-0.4.2734.tar.gz:
    https://wiki.fysik.dtu.dk/gpaw-files/gpaw-0.4.2734.tar.gz

After getting the code :ref:`create_links`.

.. _latest_development_release:

Latest development release
==========================

The latest revision can be obtained from svn::

  $ svn checkout https://svn.fysik.dtu.dk/projects/gpaw/trunk gpaw

or from the daily snapshot: `<gpaw-snapshot.tar.gz>`_.

After getting the code :ref:`create_links`.

.. note::

   The recommended checkout path is :envvar:`HOME`.

See :ref:`faq` in case of problems.

.. _create_links:

Create links
============

It is convenient to maintain several version of GPAW
with the help of links.
After downloading create the link to the requested version, e.g.:

- if retrieved from ``svn``::

   $ cd $HOME
   $ ln -s gpaw-0.9.0 gpaw

- if retrieved as ``tarball``::

   $ cd $HOME
   $ tar -xtf gpaw-0.9.0.8965.tar.gz
   $ ln -s gpaw-0.9.0.8965 gpaw

  .. note::

     The recommended installation path is :envvar:`HOME`.

When you have the code, go back to the :ref:`installationguide_developer`.
