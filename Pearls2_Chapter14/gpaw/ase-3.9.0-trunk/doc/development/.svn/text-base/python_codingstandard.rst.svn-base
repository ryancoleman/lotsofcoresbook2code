.. _python_codingstandard:

==================
Coding Conventions
==================


Importing modules
=================

In code, like the implementation of ASE, we must *not* use the
``import *`` syntax.  Import everything explicitly from exactly the
place where it's defined::

  from ase.io import read, write

We distinguish between scripts and code.  In your own scripts, it's OK
to use::

  from ase.all import *

which will give you the most used symbols.


Python Coding Conventions
=========================

Please run :ref:`pep8.py <pep8py>` and :ref:`pylint <pylint>` on your
code before committing.

The rules for the Python part are almost identical
to those used by the `Docutils project`_:

Contributed code will not be refused merely because it does not
strictly adhere to these conditions; as long as it's internally
consistent, clean, and correct, it probably will be accepted.  But
don't be surprised if the "offending" code gets fiddled over time to
conform to these conventions.

The project shall follow the generic coding conventions as
specified in the `Style Guide for Python Code`_ and `Docstring
Conventions`_ PEPs, summarized, clarified, and extended as follows:

* 4 spaces per indentation level.  No hard tabs.

* Very important:  Read the *Whitespace in Expressions and Statements*
  section of PEP8_.

* Avoid introducing `trailing whitespaces`_.

* Try to use only 7-bit ASCII, no 8-bit strings.

* No one-liner compound statements (i.e., no ``if x: return``: use two
  lines & indentation), except for degenerate class or method
  definitions (i.e., ``class X: pass`` is OK.).

* Lines should be no more than 78 characters long.

* Use "StudlyCaps" for class names.

* Use "lowercase" or "lowercase_with_underscores" for function,
  method, and variable names.  For short names, maximum two words,
  joined lowercase may be used (e.g. "tagname").  For long names with
  three or more words, or where it's hard to parse the split between
  two words, use lowercase_with_underscores (e.g.,
  "note_explicit_target", "explicit_target").  If in doubt, use
  underscores.

* Avoid lambda expressions, which are inherently difficult to
  understand.  Named functions are preferable and superior: they're
  faster (no run-time compilation), and well-chosen names serve to
  document and aid understanding.

* Avoid functional constructs (filter, map, etc.).  Use list
  comprehensions instead.

* Avoid ``from __future__ import`` constructs.  They are inappropriate
  for production code.

* Use 'single quotes' for string literals, and """triple double
  quotes""" for :term:`docstring`\ s.  Double quotes are OK for
  something like ``"don't"``.

.. _PEP8:
.. _Style Guide for Python Code: http://www.python.org/peps/pep-0008.html
.. _Docstring Conventions: http://www.python.org/peps/pep-0257.html
.. _Docutils project: http://docutils.sourceforge.net/docs/dev/policies.html#python-coding-conventions
.. _trailing whitespaces: http://www.gnu.org/software/emacs/manual/html_node/emacs/Useless-Whitespace.html

.. attention::

   Thus spake the Lord: Thou shalt indent with four spaces. No more, no less.
   Four shall be the number of spaces thou shalt indent, and the number of thy
   indenting shall be four. Eight shalt thou not indent, nor either indent thou
   two, excepting that thou then proceed to four. Tabs are right out.

                                          Georg Brandl


General advice
==============

 * Get rid of as many ``break`` and ``continue`` statements as possible.


Writing documentation in the code
=================================

Here is an example of how to write good docstrings:

  https://github.com/numpy/numpy/blob/master/doc/example.py


.. _pep8py:

Run pep8.py on your code
========================

The `pep8.py <https://github.com/jcrocholl/pep8>`_ program is
installed together with ASE :svn:`tools/pep8.py`.
It will check the PEP8_ conventions for you.  Try::

  $ pep8.py --help

.. _autopep8py:

Run autopep8.py on your code
============================

Another method of enforcing PEP8_ is using a tool such as 
`autopep8.py <https://github.com/hhatto/autopep8>`_. These tools tend to be
very effective at cleaning up code, but should be used carefully and code
should be retested after cleaning it. Try::

  $ autopep8.py --help

.. attention::

   There is a common issue with pep8 where spaces are added around the power
   operator.  Code such as "x**2" should not be changed to "x ** 2".  This
   issue is not fixed in pep8 as of the time of this writing, but a small
   `change <http://listserv.fysik.dtu.dk/pipermail/gpaw-developers/2014-October/005075.html>`_
   to autopep8 has been effective to prevent this change.

.. _pylint:

Using pylint to check your code
===============================

A pylintrc trying to follow ASE :ref:`python_codingstandard` can be found here:
:svn:`doc/development/pylintrc`


Running pylint yourself
-----------------------

Run pylint on a single file like this::

    [~]$ pylint mypythonfile.py

Run pylint on a module like this::
    
    [~]$ pylint path/to/module/root/dir

.. _epydoc:

Run epydoc on your code
=======================

Run::

  $ epydoc --docformat restructuredtext --parse-only --show-imports -v dir
