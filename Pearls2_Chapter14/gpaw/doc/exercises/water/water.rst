===========================
Basics of GPAW calculations
===========================

Atomization energy revisited
============================

The below script calculates the total energies of :mol:`H_2O`, H and O.

.. literalinclude:: h2o.py

Read the script and try to understand what it does.  A few notes:

 * In most languages it is common to declare loops with integer counters.
   In Python, a ``for`` loop always loops over elements of an
   *iterable* (in this case a list of the strings ``H2O``, ``H`` and ``O``).
   You can loop over integers, if desired, using ``for x in range(17):``.

 * The code in loops, if-statements and other code blocks is indented.  
   The loop or if-statement stops when the lines are no longer indented.
   Thus, *indentation determines control flow*.

 * In this case we conveniently load the geometry from the G2 database
   of small molecules, using the :func:`~ase.structure.molecule`
   function from ASE.

 * By setting the ``txt`` parameter, we specify a file where GPAW will save
   the calculation log.

 * The expression ``'results-%.2f.txt' % h`` inserts the value of ``h``
   in place of the *substitution code* ``%.2f`` (floating point number
   with 2 decimals).  Thus the result file name evaluates to
   ``results-0.20.txt``.  Similarly, ``'gpaw-%s-%.2f.txt' % (name, h)``
   evaluates to ``gpaw-H2O-0.20.txt`` in the first loop iteration
   (``%s`` is a substitution code for a string).

 * The call to ``open`` opens a file.  The parameter ``'w'`` signifies that
   the file is opened in write mode (deleting any previous file with that
   name!).  The calculated energies are written to this file using a print
   statement.  We could have chosen to just print the parameter as in the last
   exercise, but file handling will come in handy below.

Run the script.  You can monitor the progress by opening one of the
log files (e.g. ``gpaw.H2O.txt``).  The command :samp:`tail -f
{filename}` can be used to view the output in real-time.  The calculated
atomization energy can be found in the ``results-0.20.txt`` file.


Parallelization
===============

To speed things up, let us run the script using some more CPUs.  The
actual GPAW calculation is coded to make proper use of these extra
CPUs, while the rest of the script (everything except
``calc.get_potential_energy()``) is actually going to be run
independently by each CPU.  This matters little except when writing to
files.  Each CPU will attempt to write to the same file, probably at
the same time, producing garbled data.  We must therefore make sure
that only one process writes.  ASE provides the handy
:func:`~ase.parallel.paropen` function for just that::

  from ase.parallel import paropen
  ...
  resultfile = paropen('results-%.2f.txt' % h, 'w')

Apply the above modifications to the script and run it in parallel
e.g. on four CPUs::

  $ mpirun -np 4 gpaw-python script.py

Verify by checking the log file that GPAW is actually using multiple
CPUs.  The log file should reveal that the :mol:`H_2O` calculation
uses domain-decomposed with four domains, while the two atomic
calculations should parallelize over the two spins and two domains.


.. _convergence_checks:

Convergence checks
==================

It is essential that the calculations use a sufficiently fine grid
spacing, and that the cell is sufficiently large not to affect the
result.  For this reason, convergence with respect to these parameters
should generally be checked.  For now we shall only bother to check
the grid spacings.

Modify the above script to include a loop over different grid
spacings.  For :ref:`technical reasons <poisson_performance>`, GPAW
will always use a number of grid points divisible by four along each
direction.  Use a loop structure like::

  for symbol in [...]:
      ...
      
      for ngridpoints in [24, 28, ...]:
          h = a / ngridpoints
          calc.set(h=h)
          energy = system.get_potential_energy()
          ...


The ``set`` method can be used to change the parameters of a
calculator without creating a new one.  Make sure that the numbers of
grid points are chosen to cover :math:`0.15<h<0.25`.  While performing
this convergence check, the other parameters do not need to be
converged - you can reduce the cell size to e.g. ``a = 6.0`` to
improve performance.  You may wish to run the calculation in parallel.

 * How do the total energies converge with respect to grid spacing?
 * How does the atomization energy converge?

Total energies (maybe surprisingly) do not drop
as the grid is refined.  This would be the case in plane-wave methods,
where increase of the planewave cutoff strictly increases the quality
of the basis.  Grid-based methods rely on *finite-difference
stencils*, where the gradient in one point is calculated from the
surrounding points.  This makes the grid strictly inequivalent to a
basis, and thus not (necessarily) variational.


LCAO calculations
=================

GPAW supports an alternative calculation mode, :ref:`LCAO mode
<lcao>`, which uses linear combinations of pre-calculated atomic
orbitals to represent the wavefunctions.  LCAO calculations are much
faster for most systems, but also less precise.

Performing an LCAO calculation requires setting the ``mode`` and
normally a ``basis``::

  calc = GPAW(...,
              mode='lcao',
              basis='dzp',
              ...)

Here ``dzp`` ("double-zeta polarized") means two basis functions per
valence state plus a polarization function - a function corresponding
to the lowest unoccupied angular-momentum state on the atom.  This
will use the standard basis sets distributed with GPAW.  You can pick
out a smaller basis set using the special syntax ``basis='sz(dzp)'``.
This will pick out only one function per valence state
("single-zeta"), making the calculation even faster but less precise.

Calculate the atomization energy of :mol:`H_2O` using different basis
sets.  Instead of looping over grid spacing, use a loop over basis keywords::

  for basis in ['sz(dzp)', 'szp(dzp)', 'dzp']:
      ...
      calc = GPAW(mode='lcao',
                  basis=basis, 
                  ...)

Compare the calculated energies to those calculated in grid mode.  Do
the energies deviate a lot?  What about the atomization energy?  Is
the energy variational with respect to the quality of the basis?

LCAO calculations do not in fact produce very precise binding energies
(although these can be improved considerably by manually generating
optimized basis functions) - however the method is well suited
to calculate geometries, and for applications that require a small basis
set, such as electron transport calculations.


Plane-wave calculations
=======================

For systems with small unit-cells, it can be much faster to expand the
wave-functions in :ref:`plane-waves <manual_mode>`.  Try running a calculation
for a water molecule with a plane-wave cutoff of 350 eV using this::
    
    from gpaw import GPAW, PW
    calc = GPAW(mode=PW(350), ...)

Try to look at the text output and see if you can find the number of
plane-waves used.
