================================
Nudged elastic band calculations
================================

Self-diffusion on the Al(110) surface
-------------------------------------

.. image:: Al110slab.png
   :height: 270 px
   :alt: Al(110) surface
   :align: right

In this exercise, we will find minimum-energy paths and transition states
using the :ase:`Nudged Elastic Band <ase/neb.html>` method. Another method for
finding the transition state (i.e. the highest-energy state), the Dimer
method, will also be explored.

Take a look at the Al(110) surface shown in the picture on the right. The red
atom represents an Al adatom that can move around on the surface. The adatom
can jump along the rows (into the picture) or across the rows (to the right in
the picture).

* Which of the two jumps do you think will have the largest energy
  barrier?

The template script :download:`neb1.py` will
find the minimum-energy path for a jump along the rows. Read,
understand, and run the script.

* Make sure you understand what is going on (make a good sketch of the
  110 surface).

* View the profile of the NEB path in ase-gui. How is the shape
  (symmetric/asymmetric) and does this make sense for this process
  (when looking at the moving adatom in the simulation)?

* What is the energy barrier?

* Copy the script to ``neb2.py`` and modify it to find the barrier for
  diffusion across one of the rows.  What is the barrier for this
  process?

* Can you think of a third type of diffusion process?  Hint: It is
  called an exchange process and you can read more about it in the paper listed
  :mod:`here <ase.dimer>`.
  Find the barrier for this process, and
  compare the energy barrier with the two other ones.
  (If you give up look at :download:`neb3.py`)

* Could there be other final-image configurations for the exchange process?

.. hint::

  When opening a trajectory in :program:`ase-gui` with calculated energies, the
  default plot window shows the energy versus frame number.  To get a
  better feel of the energy barrier in an NEB calculation; choose
  :menuselection:`Tools --> NEB`. This will give a smooth curve
  of the energy as a
  function of the NEB path length, with the slope at each point
  estimated from the force.

In the NEB calculations above we knew the final states, so all we had to do
was to calculate the path between the initial state and the final state. But
in some cases we do not know the final state. Then the :mod:`Dimer method
<ase.dimer>` can be used to find the transition state. The result of a Dimer
calculation will hence not be the complete particle trajectory as in the NEB
output, but rather the configuration of the transition-state image.

The template script :download:`dimer_along.py` will find the transition-state
image of the jump along the row. Again, read, understand and run the script.

* Make sure you understand what is going on. For instance, see the trajectory
  file in ase-gui.

* Compare the transition-state images of the NEB and Dimer as viewed in ase-
  gui. Are they identical?

* What is the energy barrier? How does it compare to the one found in the NEB
  calculation?

* Do the same as above for the jump across the row and the exchange process by
  copying and modifying the Dimer script, while remembering that you have to
  give the relevant atoms a kick in a meaningful direction.
