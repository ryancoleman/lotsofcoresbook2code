.. _parallel_runs:

=============
Parallel runs
=============

.. toctree::
   :maxdepth: 1

.. _parallel_running_jobs:

Running jobs in parallel
========================

Parallel calculations are done with MPI and a special
:program:`gpaw-python` python-interpreter.

The parallelization can be done over the **k**-points, bands, spin in
spin-polarized calculations, and using real-space domain
decomposition.  The code will try to make a sensible domain
decomposition that match both the number of processors and the size of
the unit cell.  This choice can be overruled, see
:ref:`manual_parallelization_types`.

Before starting a parallel calculation, it might be useful to check how the parallelization corresponding to the given number
of processes would be done with ``--dry-run`` command line option::

  python script.py --dry-run=8

The output will contain also the "Calculator" RAM Memory estimate per process.

In order to start parallel calculation, you need to know the
command for starting parallel processes. This command might contain
also the number of processors to use and a file containing the names
of the computing nodes.  Some
examples::

  mpirun -np 4 gpaw-python script.py
  poe "gpaw-python script.py" -procs 8


Simple submit tool
==================

Instead writing a file with the line "mpirun ... gpaw-python script.py" and then submitting it to a queueing system, it is simpler to automate this::

  #!/usr/bin/env python
  from sys import argv
  import os
  options = ' '.join(argv[1:-1])
  job = argv[-1]
  dir = os.getcwd()
  f = open('script.sh', 'w')
  f.write("""\
  NP=`wc -l < $PBS_NODEFILE`
  cd %s
  mpirun -np $NP -machinefile $PBS_NODEFILE gpaw-python %s
  """ % (dir, job))
  f.close()
  os.system('qsub ' + options + ' script.sh')

Now you can do::

  $ qsub.py -l nodes=20 -m abe job.py

You will have to modify the script so that it works with your queueing
system.

.. _submit_tool_on_niflheim:

Submit tool on Niflheim
-----------------------

At CAMd, we use this submit tool: :svn:`~doc/documentation/parallel_runs/gpaw-qsub`.

Examples::

  $ gpaw-qsub -q medium -l nodes=8 -m abe fcc.py --domain-decomposition=1,2,2
  $ gpaw-qsub -q long -l nodes=6:ppn=8:xeon5570 -m abe hcp_n2.py --gpaw=blacs=1 \
    --sl_default=4,4,2 --domain-decomposition=8 --state-parallelization=2

.. tip::
  CAMd users must always remember to source the openmpi environment settings before recompiling the code. See :ref:`Niflheim`.

Alternative submit tool
=======================

Alternatively, the script gpaw-runscript can be used, try::

  $ gpaw-runscript -h

to get the architectures implemented and the available options. As an example, use::

  $ gpaw-runscript script.py 32

to write a job sumission script running script.py on 32 cpus. 
The tool tries to guess the architecture/host automatically.

By default it uses the following environment variables to write the runscript:

=============== ===================================
variable        meaning                            
=============== ===================================
HOSTNAME        name used to assing host type
PYTHONPATH      path for python
GPAW_PYTHON     where to find gpaw-python
GPAW_SETUP_PATH where to find the setups
GPAW_MAIL       where to send emails about the jobs
=============== ===================================

Writing to files
================

Be careful when writing to files in a parallel run.  Instead of ``f = open('data', 'w')``, use:

>>> from ase.parallel import paropen
>>> f = paropen('data', 'w')

Using ``paropen``, you get a real file object on the master node, and dummy objects on the slaves.  It is equivalent to this:

>>> from ase.parallel import rank
>>> if rank == 0:
...     f = open('data', 'w')
... else:
...     f = open('/dev/null', 'w')

If you *really* want all nodes to write something to files, you should make sure that the files have different names:

>>> from ase.parallel import rank
>>> f = open('data.%d' % rank, 'w')

Writing text output
===================

Text output written by the ``print`` statement is written by all nodes.
To avoid this use:

>>> from ase.parallel import parprint
>>> print 'This is written by all nodes'
>>> parprint('This is written by the master only')

which is equivalent to

>>> from ase.parallel import rank
>>> print 'This is written by all nodes'
>>> if rank == 0:
...     print 'This is written by the master only'

Note that parprint has the syntax of the print statement in 
`Python3 <http://docs.python.org/release/3.0.1/whatsnew/3.0.html>`_.

.. _different_calculations_in parallel:

Running different calculations in parallel
==========================================

A GPAW calculator object will per default distribute its work on all
available processes. If you want to use several different calculators
at the same time, however, you can specify a set of processes to be
used by each calculator. The processes are supplied to the constructor,
either by specifying an :ref:`MPI Communicator object <communicators>`,
or simply a list of ranks. Thus, you may write::

  from gpaw import GPAW
  import gpaw.mpi as mpi

  # Create a calculator using ranks 0, 3 and 4 from the mpi world communicator
  ranks = [0, 3, 4]
  comm = mpi.world.new_communicator(ranks)
  if mpi.world.rank in ranks:
      calc = GPAW(communicator=comm)
      ...

Be sure to specify different output files to each calculator,
otherwise their outputs will be mixed up.

Here is an example which calculates the atomization energy of a
nitrogen molecule using two processes:

.. literalinclude:: parallel_atomization.py

.. _manual_parallelization_types:

.. _manual_parallel:

Parallization options
=====================

In version 0.7, a new keyword called ``parallel`` was introduced to provide 
a unified way of specifying parallelization-related options. Similar to
the way we :ref:`specify convergence criteria <manual_convergence>` with the 
``convergence`` keyword, a Python dictionary is used to contain all such
options in a single keyword.

The default value corresponds to this Python dictionary::

  {'kpt':                 None,
   'domain':              None,
   'band':                1,
   'order':               'kdb',
   'stridebands':         False,
   'sl_auto':             False,
   'sl_default':          None,
   'sl_diagonalize':      None,
   'sl_inverse_cholesky': None,
   'sl_lcao':             None,
   'sl_lrtddft':          None,
   'buffer_size':         None}

In words:

* ``'kpt'`` is an integer and denotes the number of groups of k-points over which to parallelize.
  k-point parallelization is the most efficient type of parallelization for most systems
  with many electrons and/or many k-points.
  If uspecified, the calculator will choose a parallelization itself which maximizes the k-point
  parallelization unless that leads to load imbalance; in that case, it may prioritize domain
  decomposition.

* The ``'domain'`` value specifies either an integer ``n`` or a tuple
  ``(nx,ny,nz)`` of 3 integers for :ref:`domain decomposition <manual_parsize_domain>`.
  If not specified (i.e. ``None``), the calculator will try to determine the best
  domain parallelization size based on number of kpoints, spins etc.

* The ``'band'`` value specifies the number of parallelization groups to use for
  :ref:`band parallelization <manual_parsize_bands>` and defaults to one, i.e.
  no band parallelization.

* ``'order'`` specifies how different parallelization modes are nested
  within the calculator's world communicator.  Must be a permutation
  of the characters ``'kdb'`` which is the default.  The characters
  denote k-point, domain or band parallelization respectively.  The
  last mode will be assigned contiguous ranks and thus, depending on
  network layout, probably becomes more efficient.  Usually for static
  calculations the most efficient order is ``'kdb'`` whereas for TDDFT
  it is ``'kbd'``.

* The ``'stridebands'`` value only applies when band parallelization is used, and
  can be used to toggle between grouped and strided band distribution.

* If ``'sl_auto'`` is ``True``, ScaLAPACK will be enabled with automatically chosen
  parameters and using all available CPUs.

* The other ``'sl_...'`` values are for using ScaLAPACK with different parameters in different operations.
  Each can be specified as a tuple ``(m,n,mb)`` of 3 integers to
  indicate an ``m*n`` grid of CPUs and a block size of ``mb``.
  If any of the three latter keywords are not
  specified (i.e. ``None``), they default to the value of
  ``'sl_default'``. Presently, ``'sl_inverse_cholesky'`` must equal
  ``'sl_diagonalize'``.

* The ``'buffer_size'``  is specified as an integer and corresponds to
  the size of the buffer in KiB used in the 1D systolic parallel
  matrix multiply algorithm. The default value corresponds to sending all
  wavefunctions simultaneously. A reasonable value would be the size
  of the largest cache (L2 or L3) divide by the number of MPI tasks
  per CPU. Values larger than the default value are non-sensical and
  internally reset to the default value.

.. note::
   With the exception of ``'stridebands'``, these parameters all have an
   equivalent command line argument which can equally well be used to specify
   these parallelization options. Note however that the values explicitly given
   in the ``parallel`` keyword to a calculator will override those given via
   the command line. As such, the command line arguments thus merely redefine
   the default values which are used in case the ``parallel`` keyword doesn't
   specifically state otherwise.


.. _manual_parsize_domain:

Domain decomposition
--------------------

Any choice for the domain decomposition can be forced by specifying
``domain`` in the ``parallel`` keyword. It can be given in the form
``parallel={'domain': (nx,ny,nz)}`` to force the decomposition into ``nx``,
``ny``, and ``nz`` boxes in x, y, and z direction respectively. Alternatively,
one may just specify the total number of domains to decompose into, leaving
it to an internal cost-minimizer algorithm to determine the number of domains
in the x, y and z directions such that parallel efficiency is optimal. This
is achieved by giving the ``domain`` argument as ``parallel={'domain': n}``
where ``n`` is the total number of boxes.

.. tip::
   ``parallel={'domain': world.size}`` will force all parallelization to be
   carried out solely in terms of domain decomposition, and will in general
   be much more efficient than e.g. ``parallel={'domain': (1,1,world.size)}``.
   You might have to add ``from gpaw.mpi import wold`` to the script to 
   define ``world``.

There is also a command line argument ``--domain-decomposition`` which allows you
to control domain decomposition (see example at :ref:`submit_tool_on_niflheim`).

.. _manual_parsize_bands:

Band parallelization
--------------------

Parallelization over Kohn-Sham orbitals (i.e. bands) becomes favorable when
the number of bands :math:`N` is so large that :math:`\mathcal{O}(N^2)`
operations begin to dominate in terms of computational time. Linear algebra
for orthonormalization and diagonalization of the wavefunctions is the most
noticeable contributor in this regime, and therefore, band parallelization
can be used to distribute the computational load over several CPUs. This
is achieved by giving the ``band`` argument as ``parallel={'band': nbg}``
where ``nbg`` is the number of band groups to parallelize over.

.. tip::
   Whereas band parallelization in itself will reduce the amount of operations
   each CPU has to carry out to calculate e.g. the overlap matrix, the actual
   linear algebra necessary to solve such linear systems is in fact still
   done using serial LAPACK by default. It is therefor advisable to use both
   band parallelization and ScaLAPACK in conjunction to reduce this
   potential bottleneck.

There is also a command line argument ``--state-parallelization`` which allows you
to control band parallelization (see example at :ref:`submit_tool_on_niflheim`).

More information about these topics can be found here:

.. toctree::
   :maxdepth: 1

   band_parallelization/band_parallelization

.. _manual_ScaLAPACK:

ScaLAPACK
--------------------
Calculations with nbands > 500 will benefit from ScaLAPACK (otherwise
the default serial LAPACK should be used). The ScaLAPACK parameters
are defined either using the aformentioned ``'sl_...'`` entry in the parallel
keyword dictionary or using a command line argument,
e.g. ``--sl_default=m,n,mb``.  A reasonbly good guess for these 
parameters on most systems is related to the numbers of bands. 
We recommend::

  mb = 64
  m = floor(sqrt(nbands/mb))
  n = m

There are a total of four ``'sl_...'`` keywords. Most people will be
fine just using ``'sl_default'``. Here we use the same 
ScaLAPACK parameters in three different places: i) general eigensolve
in the LCAO intilization ii) standard eigensolve in the FD calculation and 
iii) Cholesky decomposition in the FD calculation. It is currently
possible to use different ScaLAPACK parameters in the LCAO
intialization and the FD calculation by using two of the ScaLAPACK 
keywords in tandem, e.g::

   --sl_lcao=p,q,pb --sl_default=m,n,mb

where ``p``, ``q``, ``pb``, ``m``, ``n``, and ``mb`` all 
have different values. The most general case is the combination
of three ScaLAPACK keywords, e.g::

   --sl_lcao=p,q,pb --sl_diagonalize=m,n,mb  --sl_inverse_cholesky=r,s,rb

however, we do not presently support ``m != r``, ``n != s``,  and
``mb != rb``.  We may implement this in the future.

