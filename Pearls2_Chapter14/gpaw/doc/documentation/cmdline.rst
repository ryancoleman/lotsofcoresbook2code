.. program:: gpaw

.. index:: gpaw, command line tool

.. _command line tool:

======================
GPAW command line tool
======================

.. highlight:: bash

The :program:`gpaw` command line tool is nothing more than a shortcut
to the :program:`ase` tool, with the calculator explicitely set to GPAW.
Read the ASE documentation here: `ASE command line tool`_.  Instead of
specifying the GPAW calculator like here::

    $ ase gpaw bulk Li -k 2,2,2

you can use the :program:`gpaw` command directly::

    $ gpaw bulk Li -k 2,2,2


.. _ASE command line tool: https://wiki.fysik.dtu.dk/ase/ase/cmdline.html


Special GPAW options
====================

GPAW specific options:

--parameter-file=FILE
                    Read GPAW parameters from file.
-S, --show-text-output
                    Send text output from calculation to standard out.
-W MODE, --write-gpw-file=MODE
                    Write gpw file.
