.. _summerschool08:

=======================
CAMd Summer school 2008
=======================

Databar
=======

Setting up your UNIX environment
--------------------------------

The first time you use the databar computers, you must do this:

.. highlight:: bash

::

  $ ~jjmo/summerschool/setup.sh
  $ source ~/.bashrc

That will set up the environment for you so that you can use ASE, GPAW, VMD and matlpotlib.

**Warning** runnig :command:`~jjmo/summerschool/setup.sh` owervrites
users ~/.bashrc ~/.emacs ~/.pythonrc ~/.vmdrc and ~/.matplotlib directory.

Notes
-----

* Useful links: Userguides_ FAQ_ Unix_ USB-sticks_

* Editors: emacs, vim, nedit (MS Windows/Macintosh-like environment). Python syntax

* Printer: gps1-308. Terminal: lp -d gps1-308 filename

* E-mail client:
  Thunderbird is the default mail client in the databar and configured  
  with your summer school e-mail (camd0??@student.dtu.dk).

* To open a pdf-file: acroread filename

Niflheim
========

Frontend nodes
--------------

Log in to niflheim::

  ssh school1.fysik.dtu.dk

or::

  ssh school2.fysik.dtu.dk.

Same password as handed out for the databar. Please use school1 if the
number in your userid is odd and school2 if it is even.

Copying files from gbar to niflheim
-----------------------------------

You can copy files from the Gbar to niflheim with ``scp``. If you are on 
niflheim::

    scp hald.gbar.dtu.dk:path/filename .

will copy ``filename`` to your present location. Likewise::

    scp school1.fysik.dtu.dk:path/filename .

will copy ``filename`` from Niflheim to your present location at the Gbar.

GPAW
----

Use the :command:`gpaw-qsub.py` command to submit GPAW jobs to the queue.


SIESTA
------

Siesta is installed on Niflheim, so you need to log in to the Niflheim
front-end nodes as described above in the Niflheim section.
Furthermore you have to set two environment variables by adding the
following two lines to your ~/.bashrc file::

  export SIESTA_PP_PATH=~mvanin/asesiesta
  export SIESTA_SCRIPT=~mvanin/asesiesta/run_siesta.py  

and source it by typing::

  $ source ~/.bashrc

To submit a job to Niflheim, use the ``qsub`` command::

  $ qsub -l nodes=1:ppn=1:switch5 filename.py


Octopus
-------

Octopus_ is installed on the 'q' opteron nodes on Niflheim. The way to
run jobs is the following: Create inp file in the working directory as
described in the tutorial_, and then run
:svn:`~doc/exercises/octopus_run.py`. To use various octopus utilities such
as ``oct-cross-section`` and ``oct-broad`` you need to do::

  source /usr/local/openmpi-1.2.5-pathf90/bin/mpivars-1.2.5.sh

first. Submitting jobs to the queue is done by::

  qsub -l nodes=2:ppn=4:switch5 octopus_run.py


.. _Userguides: http://www.gbar.dtu.dk/index.php/Category:User_Guides
.. _FAQ: http://www.gbar.dtu.dk/index.php/General_use_FAQ
.. _Unix: http://www.gbar.dtu.dk/index.php/UNIX
.. _USB-sticks: http://www.gbar.dtu.dk/index.php/USBsticks
.. _Octopus: http://www.tddft.org/programs/octopus/wiki/index.php/
.. _tutorial: http://www.tddft.org/programs/octopus/wiki/index.php/Tutorial

