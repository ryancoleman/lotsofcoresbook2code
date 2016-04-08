.. _bugs:

=======================
Bugs in the latest GPAW
=======================

Bugs and feature proposals are managed using `Trac <https://trac.fysik.dtu.dk/projects/gpaw/>`_. If you have
a problem with GPAW, and you are already subscribed to the :ref:`mailing_lists`, this is the place to handle:

* Bugs in the code
* Enhancement proposals
* Problems with the documentation

Remember to check the list of `tickets <https://trac.fysik.dtu.dk/projects/gpaw/report/1>`_ for known issues first.

------------------
Handling segfaults
------------------

Segmentation faults are probably the hardest type of runtime error to track down, but they are also quite
common during the *unstable* part of the release cycle. As a rule of thumb, if you get a segfault, start
by checking that all array arguments passed from Python to C functions have the correct shapes and types.

Apart from appending ``--debug`` to the command line arguments when running ``python`` or ``gpaw-python``,
please familiarize yourself with the :ref:`debugging tools <debugging>` for the Python and C code.

If you experience segfaults or unexplained MPI crashes when running GPAW in parallel, it is recommended to
try a :ref:`custom installation <install_custom_installation>` with a debugging flag in ``customize.py``::

  define_macros += [("GPAW_MPI_DEBUG",1)]


----------------------
Common sources of bugs
----------------------

* General:

  - Look for XXX, TODO, !!! or ??? in the source code.
  - Code which is not used very often will tend to stop working.
  - Be weary of copy/paste errors. Avoid code duplication if you can.
  - Numerics default type on an alpha is Int64. Use long instead of int.
  - Elements of NumPy arrays are C ordered, BLAS and LAPACK routines expect Fortran ordering.  

.. spacer

* Python:

  - Sometimes you need a copy and not a reference.

  - Always give contiguous arrays to C functions. If ``x`` is contiguous with
    ``dtype=complex``, then ``x.real`` is non-contiguous of ``dtype=float``.

  - Giving array arguments to a function is a *carte blanche* to alter the data::

      def double(a):
          a *= 2
          return a

      x = np.ones(5)
      print double(x) # x[:] is now 2.

    NumPy arrays can be made read-only to prevent accidental changes::

      x.flags.writable = False
      x *= 2

  - Forgetting a ``n += 1`` statement in a for loop::

      n = 0
      for thing in things:
	  thing.do_stuff(n)
	  n += 1

    Use this instead::

      for n, thing in enumerate(things):
	  thing.do_stuff(n)

  - Indentation errors like this one::

     if ok:
         x = 1.0
     else:
         x = 0.5
         do_stuff(x)

    where ``do_stuff(x)`` should have been reached in both cases.
    Emacs: always use ``C-c >`` and ``C-c <`` for shifting in and out
    blocks of code (mark the block first).

  - Don't use mutables as default values::

     class A:
         def __init__(self, a=[]):
             self.a = a # All instances get the same list!

  - There are subtle differences between ``x == y`` and ``x is y``.

  - If ``H`` is a numeric array, then ``H - x`` will subtract ``x``
    from *all* elements - not only the diagonal, as in Matlab!


* C:

  - Try building GPAW from scratch.
  - Typos like ``if (x = 0)`` which should have been ``if (x == 0)``.
  - Remember ``break`` in switch-case statements.
  - Check ``malloc-free`` pairs. Test for :ref:`memory leaks <memory_leaks>` by repeating the call many times.
  - Remember to update reference counts of Python objects.
  - *Never* put function calls inside ``assert``'s.  Compiling with
    ``-DNDEBUG`` will remove the call.

