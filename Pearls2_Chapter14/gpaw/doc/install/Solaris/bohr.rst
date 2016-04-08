.. _bohr:

================
bohr.gbar.dtu.dk
================

Here you find information about the the system
`<http://www.gbar.dtu.dk/index.php/Hardware>`_.

Follow instructions from `<http://www.gbar.dtu.dk/index.php/GridEngine>`_ to create :file:`~/.grouprc`.

Use ``python2.5`` and enable ``numpy``/``matplotlib``::

 alias python="LD_PRELOAD=/opt/csw/lib/libncurses.so /usr/local/gbar/cswbin/python2.5"
 export PYTHONPATH=/usr/local/gbar/lib/pythonmodules

To build gpaw add to the ``gpaw/customize.py``:

  library_dirs += ['/opt/csw/lib']

Download `MPIscript.sh <http://www.hpc.dtu.dk/GridEngine/MPIscript.sh>`_ and edit it, so it resembles::

 #!/bin/sh 
 # (c) 2000 Sun Microsystems, Inc.
 # ---------------------------
 # General options
 #
 #$ -S /bin/sh
 #$ -o $JOB_NAME.$JOB_ID.out
 #$ -e $JOB_NAME.$JOB_ID.err
 # -M User@Domain
 # -m es
 # ---------------------------
 # Execute the job from  the  current  working  directory
 #$ -cwd
 #
 # Parallel environment request
 # ---------------------------
 # do not change the following line
 #$ -l cre
 #
 #      PE_name  CPU_Numbers_requested
 ##$ -pe HPC      2
 # ------------------------------- Program_name_and_options
 LD_LIBRARY_PATH=/opt/csw/lib:${LD_LIBRARY_PATH}
 export LD_LIBRARY_PATH
 gpawpython=~mdul/gpaw/build/bin.solaris-2.10-sun4u-2.5/gpaw-python
 /appl/hgrid/current/bin/mprun -np $NSLOTS $gpawpython gpaw/examples/H.py
 # ---------------------------

Submit jobs like this::

  qsub -N test -pe HPC 2 MPIscript.sh


.. note::

   All scripts making use of ``#!/usr/bin/env python`` must be changed
   to use ``#!/usr/bin/env python2.5`` instead.

.. _bohr_gbar_dtu_dk:
