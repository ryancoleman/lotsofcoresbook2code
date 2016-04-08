==================
Command line tools
==================

This proposal tries to clean up the command line tools that we
distribute with ASE.


Current status
==============

ASE currently has 7 command line tools that we install in
:file:`/usr/bin`:

=========================  ===============================================
command                    description
=========================  ===============================================
:command:`ag`              ASE's GUI
:command:`ASE2ase`         Translate old ASE-2 code to ASE-3 style
:command:`testase`         Run tests
:command:`ase`             First attempt to create a command that can do
                           *everything*
:command:`asec`            Second attempt
:command:`foldtrajectory`  Wrap atoms outside simulation box to inside box
:command:`trajectoryinfo`  Write information about trajectory file
=========================  ===============================================


Proposed set of command line tools
==================================

In the future things will look like this:

====================  =======================================
command               description
====================  =======================================
:command:`ase-gui`    ASE's GUI
:command:`ase-info`   Write information about files
:command:`ase-test`   Run tests
:command:`ase-build`  Build simple molecule or bulk structure
:command:`ase-run`    Run calculations with ASE's calculators
:command:`ase-db`     Put stuff into or query database
====================  =======================================


Comments
========

:command:`ag`:

    Renamed to :command:`ase-gui`.

:command:`ASE2ase`:

    Removed --- no longer needed.

:command:`testase`:

    Renamed to :command:`ase-test`.  Alternative::

        python setup.py test

:command:`ase` and :command:`asec`:

    Replaced by new commands :command:`ase-build` and
    :command:`ase-run`.  The old :command:`ase` command is documented
    :ref:`here <command line tool>` and is hopefully not used a lot
    since we propose to get rid of it.

:command:`foldtrajectory`:

    Too specialized to deserve its own command.  Use::

        python -m ase.md.foldtrajectory

    instead.

:command:`trajectoryinfo`:

    Replaced by new more general command :command:`ase-info` that can
    pull out information from anything that ASE can read.


Naming convention
=================

Any suggestions for better names or are the proposed ones OK?  The
good thing about using :command:`ase-something` for all is that it is
consistent and if you know one command, you will maybe discover the
other ones when you do tab-completion.


Implementation details
======================

* Should we use the very nice :mod:`argparse` module, which is the
  future but only available in Python 2.7, or should we stick with the
  old and deprecated :mod:`optparse` module?
