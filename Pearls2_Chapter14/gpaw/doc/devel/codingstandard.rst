.. _codingstandard:

==================
Coding Conventions
==================

Python Coding Conventions
=========================

Follow :ase:`ASE's rules <development/python_codingstandard.html>`.

C-code
======

Code C in the C99 style::

  for (int i = 0; i < 3; i++) {
      double f = 0.5;
      a[i] = 0.0;
      b[i + 1] = f * i;
  }

and try to follow PEP7_.

Use **M-x c++-mode** in emacs.

.. _PEP7: http://www.python.org/dev/peps/pep-0007
