=================
Density of states
=================

Take a look at the :svn:`~doc/exercises/dos/dos.py` program and try to
get a rough idea of what it can do for you.  Use it to plot the
density of states (DOS) for the three Fe configurations from the
:ref:`iron_exercise` exercise (on the *x*-axis you have the energy
relative to the Fermilevel).

* Do the DOS plots integrate to the correct numbers? (i.e.
  number of bands).

* The DOS for the anti-ferromagnetic phase looks a bit like that for
  the non-magnetic phase - is it magnetic at all?!  Calculate
  the magnetization like this:

  .. literalinclude:: magnetization.py
    
  and :ref:`look at it <iso>`.

* Calculate the DOS for bulk Aluminum and compare it
  (qualitatively) to the DOS for the non-magnetic calculation. The DOS
  for a simple metal has this shape: *g*\ (*E*) ~ *E*\ :sup:`1/2`.  Explain
  the qualitative difference.

* Plot also the DOS for bulk Si and the CO molecule.  Identify the
  bandgap between valence and conduction bands for Si and the
  HOMO-LUMO gap for CO. Make sure that your **k**-point mesh for
  Si is dense enough to sample the band structure.


Projected Density of states (PDOS)
----------------------------------

The projected density of states is useful for for analyzing chemical
bonding. There exist several studies where the density projected onto
the d states of a given surface atom is used. This short excercise
demonstrates how to construct the PDOS of Fe.

We will get a feel for the local density of states by plotting the
PDOS for the ferro-magnetic Fe crystal.  Look at
:svn:`~doc/exercises/dos/pdos.py`. Use it to plot the s, p,
and d-states on one of the Fe atoms.
