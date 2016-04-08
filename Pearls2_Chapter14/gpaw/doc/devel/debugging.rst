.. _debugging:

=========
Debugging
=========

Python debugging
================

Even though some debugging can done just with print statements, a real debugger offers several advantages. It is possible, for example, to set breakpoints in certain files or functions, execute the code step by step, examine and change values of variables. Python contains a standard debugger *pdb*. A script can be started under the debugger control as *python -m pdb script.py* (python 2.4) or *python /path_to_pdb/pdb.py script.py* (python 2.3). Now, before the execution of the script starts one enters the debugger prompt. The most important debugger commands are:

h(elp) [command]

b(reak) [[filename:]lineno|function[, condition]]

  Set a breakpoint.

s(tep)

  Execute the current line, stop at the first possible occasion (either in a function that is called or on the next line 
  in the current function).

n(ext)

  Continue execution until the next line in the current function is reached or it returns. 

r(eturn)

  Continue execution until the current function returns.

c(ont(inue))

  Continue execution, only stop when a breakpoint is encountered. 

l(ist) [first[, last]]

  List source code for the current file. 

p expression

  Evaluate the expression in the current context and print its value. Note: "print" can also be used, but is not a   
  debugger command -- this executes the Python print statement

Most commands can be invoked with only the first letter. A full list of 
all the commands and their explanation can be found in the `Python debugger (PDB)
documentation <http://docs.python.org/lib/module-pdb.html>`_.


An example session might look like::

  corona1 ~/gpaw/trunk/test> python -m pdb H.py
  > /home/csc/jenkovaa/gpaw/trunk/test/H.py(1)?()
  -> from gpaw import GPAW
  (Pdb) l 11,5
   11     hydrogen.SetCalculator(calc)
   12     e1 = hydrogen.GetPotentialEnergy()
   13
   14     calc.Set(kpts=(1, 1, 1))
   15     e2 = hydrogen.GetPotentialEnergy()
   16     equal(e1, e2)
  (Pdb) break 12
  Breakpoint 1 at /home/csc/jenkovaa/gpaw/trunk/test/H.py:12
  (Pdb) c

    ... output from the script...

  > /home/csc/jenkovaa/gpaw/trunk/test/H.py(12)?()
  -> e1 = hydrogen.GetPotentialEnergy()
  (Pdb) s
  --Call--
  > /v/solaris9/appl/chem/CamposASE/ASE/ListOfAtoms.py(224)GetPotentialEnergy()
  -> def GetPotentialEnergy(self):
  (Pdb) p self
  [Atom('H', (2.0, 2.0, 2.0))]


Emacs has a special mode for python debugging which can be invoked as *M-x pdb*. After that one has to give the command to start the debugger (e.g. python -m pdb script.py). Emacs opens two windows, one for the debugger command prompt and one which shows the source code and the current point of execution. Breakpoints can be set also on the source-code window.

C debugging
===========

First of all, the C-extension should be compiled with the *-g* flag in order to get the debug information into the library. 
Also, the optimizations should be switched off which could be done in :ref:`customize.py <install_custom_installation>` as::

   extra_link_args += ['-g']
   extra_compile_args += ['-O0', '-g']

There are several debuggers available, the following example session applies to *gdb*::

  sepeli ~/gpaw/trunk/test> gdb python
  GNU gdb Red Hat Linux (6.1post-1.20040607.52rh)
  (gdb) break Operator_apply
  Function "Operator_apply" not defined.
  Make breakpoint pending on future shared library load? (y or [n]) y

  Breakpoint 1 (Operator_apply) pending.
  (gdb) run H.py
  Starting program: /usr/bin/python2.4 H.py

    ... output ...
  
  Breakpoint 2, Operator_apply (self=0x2a98f8f670, args=0x2a9af73b78)
    at c/operators.c:83
  (gdb)

One can also do combined C and python debugging by starting the input
script as ``run -m pdb H.py`` i.e::

  sepeli ~/gpaw/trunk/test> gdb python
  GNU gdb Red Hat Linux (6.1post-1.20040607.52rh)
  (gdb) break Operator_apply
  Function "Operator_apply" not defined.
  Make breakpoint pending on future shared library load? (y or [n]) y

  Breakpoint 1 (Operator_apply) pending.
  (gdb) run -m pdb H.py
  Starting program: /usr/bin/python2.4 -m pdb H.py
  [Thread debugging using libthread_db enabled]
  [New Thread -1208371520 (LWP 1575)]
  > /home/jenkovaa/test/H.py(1)?()
  -> from gpaw import GPAW
  (Pdb)


The basic gdb commands are the same as in pdb (or vice versa). Full documentation
can be found in the `GDB user manual <http://www.gnu.org/software/gdb/documentation/>`_.
Apart from the commands mentioned earlier, a few are worthy of mention here:

backtrace [n | full]

   Print a backtrace of the entire stack: one line per frame for all frames in the stack
   ``full`` prints the values of the local variables also. ``n`` specifies the number
   of frames to print

jump linespec

   Resume execution at line ``linespec`` i.e. at the given location in the
   corresponding source code. Any location of the type ``filename:linenum``
   will do, but the results may be bizarre if ``linespec`` is in a different
   function from the one currently executing.

tbreak [[filename:]lineno|function[, condition]]

   Set a breakpoint similar to how ``break`` operates, but this type of breakpoint
   is automatically deleted after the first time your program stops there.

p(rint) expr

   Inquire about the symbols (names of variables, functions and types) defined
   in a compiled program. ``expr`` may include calls to functions in the program
   being debugged. Can also be used to evaluate more complicated expressions
   or referring to static variables in other source files as ``'foo.c'::x``.   


.. hint::

   Emacs can be used also with gdb. Start with *M-x gdb* and then continue
   as when starting from the command line.

.. _memory_leaks:

Tracking memory leaks
---------------------

Although a C-extensions runs fine, or so it seems, reference counting of Python
objects and matching calls to ``malloc`` and ``free`` may not always be up to par.
Frequently, the symptom of such disproportions is all too clear, resulting in
segmentation faults (i.e. ``SIGSEGV``) e.g. when a memory address is accessed
before it has been allocated or after is has been deallocated. Such situations
can be debugged using *gdb* as described above.

.. note::

   Please refer to the Python/C API Reference Manual or the unofficial (but helpful)
   introduction to `reference counting in Python <http://edcjones.tripod.com/refcount.html>`_.

On the other hand, neglecting the deallocation or forgetting to decrease the
reference count of a Python object will lead to a build-up of unreachable
memory blocks - a process known as memory leakage. Despite being non-critical
bugs, severe memory leaks in C-code will eventually bring all computations to
a halt when the program runs out of available memory.

Suppose you have written a Python script called ``test.py`` which appears to
suffer from memory leaks. Having build GPAW with the *-g* flag as described,
tracking down the source of the memory leak (in this case line 123 of ``myfile.c``)
can be done using Valgrind_ as follows::

   sepeli ~/gpaw/trunk/test> valgrind --tool=memcheck --leak-check=yes \
   --show-reachable=yes --num-callers=20 --track-fds=yes gpaw-python test.py

   ==16442== 6,587,460 bytes in 29,943 blocks are definitely lost in loss record 85 of 85
   ==16442==    at 0x40053C0: malloc (vg_replace_malloc.c:149)
   ==16442==    by 0x5322831: ???
   ==16442==    by 0x8087BD5: my_leaky_function (myfile.c:123)

Note that Valgrind_ is more than just a memory profiler for C; it provides an
entire instrumentation framework for building dynamic analysis tools and thus
includes other debugging tools, e.g. a heap/stack/global array overrun detector.

.. _Valgrind: http://valgrind.org


.. _parallel_debugging:

Parallel debugging
==================

Debugging programs that are run in parallel with MPI is not as straight forward
as in serial, but many of the same tools can be used (e.g. GDB and Valgrind).
Note that one cannot use the Python debugger as described above because GPAW
requires that a custom Python interpreter is built with the necessary MPI bindings.

There are probably numerous ways to debug an MPI application with GDB, and experimentation
is strongly encouraged, but the following method is recommended for interactive debugging.
This approach builds upon advice in Open MPI's FAQ `Debugging applications in parallel
<http://www.open-mpi.org/faq/?category=debugging#serial-debuggers>`_, but is adapted for use
with Python on a GNU/Linux development platform. Prepend the following to your script::

   import os, sys, time, math
   from gpaw.mpi import world
   from gpaw import get_gpaw_python_path
   gpaw_python_path = os.path.join(get_gpaw_python_path(), 'gpaw-python')
   ndigits = 1 + int(math.log10(world.size))
   assert os.system('screen -S gdb.%0*d -dm gdb %s %d' \
       % (ndigits, world.rank, gpaw_python_path, os.getpid())) == 0
   time.sleep(ndigits)
   world.barrier()

This runs ``gdb /path/to/gpaw-python pid`` from within each instance of the custom Python
interpreter and detaches it into a `screen <http://www.gnu.org/software/screen/>`_ session
called ``gdb.0`` for rank 0 etc. You may now resume control of the debugger instances by
running ``screen -rd gdb.0``, entering `c` to continue and so forth for all instances.

.. hint::
   Run ``screen -ls`` to get an overview of running sessions.
   Enable logging of an attached session with Ctrl+a H (capital H).
   Use Ctrl+a Ctrl+d to detach a session but leave it running.

.. note::
   This approach only works if the problem you're trying to address occurs *after* the
   GPAW executable has been loaded. In the alternate case, it is recommended to debug
   a single instance of the parallel program with the usual serial methods first.

For details on using Valgrind on parallel programs, please refer to the online manual
`Debugging MPI Parallel Programs with Valgrind <http://valgrind.org/docs/manual/mc-manual.html#mc-manual.mpiwrap>`_
