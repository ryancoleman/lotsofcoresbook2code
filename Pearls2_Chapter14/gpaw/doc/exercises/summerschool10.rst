.. _summerschool10:

=======================
CAMd Summer school 2010
=======================

When you log into the databars, you can select various desktop
environments.  We recommend CDE or IceWM, but any will do.

Exercises will make use of terminals.  Choose :menuselection:`Terminal
--> Terminal (mrxvt)` from the GBar menu, accessible by
right-clicking on the desktop in CDE/IceWM.

Setting up your UNIX environment
--------------------------------

The first time you use the databar computers, you must configure your
environment.  Run the commands:

.. highlight:: bash

::

  $ mv ~/.bashrc ~/old.bashrc
  $ echo source ~ashj/summerschool-env/gbar-gpaw.rc > ~/.bashrc
  $ source ~/.bashrc

This will set up the environment for you so that you can use ASE,
GPAW, VMD and matplotlib.

Running GPAW calculations
-------------------------

GPAW calculations are written as Python scripts, which can be run with
the command::

  $ python filename.py

If the calculation lasts more than a few seconds, submit it to the
queue instead of running it directly::

  $ gpaw-qsub filename.py

This will allow the script to be executed on a different host, so the
jobs will be distributed efficiently even if many users logged on to
the same computer.  You can run jobs in parallel, using more CPUs for
increased speed, by specifying e.g. 4 CPUs like this::

  $ gpaw-qsub -pe 4 filename.py

The ``qstat`` or :samp:`qstat -u {USERNAME}` commands can be used to
monitor running jobs, and :samp:`qdel {JOB_ID}` to delete jobs if
necessary.


Notes
-----

* Editor: Several editors are available including emacs, vim and gedit.

* Printer: ``gps1-308``. Terminal: :samp:`lp -d gps1-308 {filename}`.  The
  printer is located in databar 15, the middle of the three databars.

* To open a pdf-file: :samp:`evince {filename}`

* How to `use USB sticks <http://www.gbar.dtu.dk/wiki/USB_Access>`_.

* The normal tilde (~) key combination is not functional on the
  databar computers.  Use :kbd:`Alt Graph + 5` to type a tilde.
