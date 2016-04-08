.. _summerschool14:

=======================
CAMd Summer School 2014
=======================


Logging in to the databar terminals
===================================

When you log into the databars, you can select various desktop
environments. Please choose IceWM, the other choices will cause trouble!

Exercises will make use of terminals.  It is important that you open
the terminal on one of the Linux servers, not on the Sun server
running the graphical user interface, as GPAW will not work on it.
Choose :menuselection:`Terminal --> xterm on Linux` from the GBar
menu, accessible by right-clicking on the desktop in IceWM.  Or press
the Terminal icon on the bottom toolbar.

Setting up your UNIX environment
================================

The first time you use the databar computers, you must configure your
environment.  Run the commands:

.. highlight:: bash

::

  $ cd
  $ mv .bashrc old.bashrc
  $ cp ~jmgla/camd.bashrc .bashrc
  $ source ~/.bashrc

This will set up the environment for you so that you can use ASE,
GPAW, VMD and matplotlib.

Note that the filename ``.bashrc`` starts with a period, making it a
hidden file in Unix.  Also note that the tilde (~) key on the databar
computers is in an unusual place: :kbd:`Alt Graph  5` 

Running GPAW calculations
=========================

**Warning** do not use spaces in the directory/file names!

GPAW calculations are written as Python scripts, which can be run with
the command::

  $ python filename.py

If the calculation lasts more than a few seconds, submit it to the
queue instead of running it directly::

  $ gpaw-qsub filename.py

This will allow the script to be executed on a different host, so the
jobs will be distributed efficiently even if many users logged on to
the same computer.  The output from your script will be written to the
files :samp:`filename.py.oNNNNNN` and :samp:`filename.py.eNNNNNN`
where :samp:`NNNNNN` is a job number.  Normal output (stdout) goes to
the :samp:`.oNNNNNN` file, whereas error messages (stderr) goes to
:samp:`.eNNNNNN`.  Unlike some queue systems these files appear when
the job starts, so you can follow job progress by looking at them.

You can run jobs in parallel, using more CPUs for
increased speed, by specifying e.g. 4 CPUs like this::

  $ gpaw-qsub -p 4 filename.py

The ``qstat`` or :samp:`qstat -u {USERNAME}` commands can be used to
monitor running jobs, and :samp:`qdel {JOB_ID}` to delete jobs if
necessary.  On the joblist from ``qstat``, you can find the JOB_ID.
You can also see the status of the jobs, Q means queued, R means
running, C means completed (jobs remain on the list for a while after
completing).


Notes and hints
===============

* Editor: Several editors are available including emacs, vim and gedit.

* Printer: There is a printer in each databar, the name is written on
  the printer. To use it from a terminal: :samp:`lp -d {printename}
  {filename}`.  Read more about printing `here <http://www.gbar.dtu.dk/wiki/Printing>`_.

* To open a pdf-file: :samp:`evince {filename}`

* How to `use USB sticks <http://www.gbar.dtu.dk/wiki/USB_Access>`_.

* The normal tilde (~) key combination is not functional on the
  databar computers.  Use :kbd:`Alt Graph + 5` to type a tilde.


Using your own laptop
=====================

If you wish to use your own laptop to log into the databar, that is
indeed possible.  You can either log in via SSH (secure shell) or
using the ThinLinc client.  Note that most likely **we cannot help you
getting it to work** if something on your laptop causes trouble.


Linux and Mac laptops
---------------------

You need to open a terminal and log in to ``login.gbar.dtu.dk``.  From
there, you log onto one of the Linux hosts.

.. highlight:: bash

::

  $ ssh -X USERNAME@login.gbar.dtu.dk
    ( ... message of the day is printed ... )
  $ linuxsh -X

Windows laptops
---------------

As there is no X server running on a windows laptop, you either have
to install one, or use a full-screen client such as ThinLinc.  

You can download ThinLinc from `Cendio's webpage`_.  Information about
how to connect with ThinLinc is available from the `DTU databar wiki`_.

Note that we recommend going into Options, Screen and disable full
screen mode.  Either set the resolution to *Near current screen
size*, or if that still gives a too large window, set the size
manually.  When logging in with thinlinc, you need to log in to the
host ``thinlinc.gbar.dtu.dk`` and proceed as if you logged into a
databar terminal.

An more pleasant (but more complicated) alternative to ThinLinc is to
install an X server.  Help on doing this can be found on the `Niflheim wiki`_.

When loggin in using an X server, you should log in to
``login.gbar.dtu.dk`` and proceed to the Linux hosts with the
command::

  $ linuxsh -X


.. _`Cendio's webpage`: http://www.cendio.com/downloads/clients/
.. _`DTU databar wiki`: http://www.gbar.dtu.dk/wiki/Thinlinc
.. _`Niflheim wiki`: https://wiki.fysik.dtu.dk/niflheim/X11_on_Windows


