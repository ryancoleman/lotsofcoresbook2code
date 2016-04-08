.. _devel:

===============
ASE development
===============

As a developer, you should subscribe to all ASE related
:ref:`mailing_lists`.


Development topics
==================

.. toctree::

   releasenotes
   contribute
   versioncontrol
   python_codingstandard
   writing_documentation_ase
   calculators
   making_movies
   newrelease
   tests
   buildbot
   translate
   todo
   py3k
   proposals/proposals


.. _devel_passwd_create:

Creating an encrypted password for SVN access
=============================================

Use this command::

  htpasswd -nm <your-desired-user-name>

and type a `good password <http://xkcd.com/936/>`__ twice.  The
encrypted password will be printed on the screen.

Alternatively, you can use::

  openssl passwd -apr1

followed by your password (twice).
