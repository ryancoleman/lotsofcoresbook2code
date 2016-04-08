.. _diffusion_exercise:

=========================================
Diffusion of gold atom on Al(100) surface
=========================================

In this ASE-tutorial:

* :ase:`Diffusion of gold atom on Al(100) surface
  <tutorials/constraints/diffusion.html>`


the energy barrier for diffusion of a gold atom on an Al(100) surface
is calculated using a semi-empirical EMT potential.  In this
exercise, we will try to use DFT and GPAW.

* Run the script from the ASE-tutorial above and use the graphical
  representation to get good initial guesses for the height of the
  gold atom in the initial and transition states (hollow and bridge
  sites).

The PAW setups for both Al and Au are quite smooth, so we can try with a low
plane-wave cutoff of 200 eV.  For a quick'n'dirty
calculation we can do with just a :math:`2 \times 2` sampling of the
surface Brillouin zone.  Use these parameters for the DFT
calculations::

  calc = GPAW(mode=PW(200), kpts=(2, 2, 1), xc='PBE')

In order to speed up the calculation, use only two frozen Al(100) layers.

* Calculate the energy of the initial and final states.  Start from
  this script: :download:`initial.py`.  Do we need to apply any
  constraint to the gold atom?

* What is the PBE energy barrier? (Do not repeat the ASE-tutorial with
  GPAW, but simply relax the gold atom at the transition state and use
  the total energy differences)

* Can both initial and transition state calculations be done with only
  one **k**-point in the irreducible part of the Brillouin zone?

* Try to repeat the EMT calculations with two frozen Al(100) layers.


Making Python Tool Boxes
========================

A science project (like the one you are going to make), will often
contain some repeated and similar sub tasks like loops over different
kind of atoms, structures, parameters etc.  As an alternative to a
plethora of similar Python scripts, made by *copy+paste*, it is
advantageous to put the repeated code into tool boxes.

Python supports such tool boxes (in Python called modules): put any
Python code into a file :file:`stuff.py` then it may be used as a tool box
in other scripts, using the Python command: ``from stuff import
thing``, where ``thing`` can be almost anything.  When Python sees
this line, it runs the file :file:`stuff.py` (only the first time) and
makes ``thing`` available.  Lets try an example:

* In file :file:`stuff.py`, put::

    constant = 17
    def function(x):
        return x - 5

* and in file :file:`program.py`, put::

    from stuff import constant, function
    print 'result =', function(constant)

* Now run the script :file:`program.py` and watch the output.

You can think of ASE and GPAW as big collections of modules, that we
use in our scripts.


Writing an adsorption script
============================

As a non-trivial example of a Python module, try to write a function:

.. function:: aual100(site, height)

The *site* argument should be one of the strings that the
:func:`ase.lattice.surface.fcc100` function accepts: ``'ontop'``,
``'hollow'`` or ``'bridge'``.  The *height* argument is the height above the
Al layer.  The function must relax a  gold atom at *site*, return the energy
and write ``<site>.txt``,  ``<site>.traj``, and ``<site>.gpw`` files. Start
from  :download:`initial.py` and make the  relevant changes.

* You could have used this functions to calculate the energy barrier
  above.  Use it to calculate the energy in the ontop site::

    e_ontop = aual100('ontop', 2.2)

* What seems to determine the relative energetic ordering of the three sites?

* Suppose now that an Au atom diffuses from one hollow to a
  neighboring hollow site at the surface.  Assuming a prefactor of 10\
  :sup:`13`/sec, how often does the diffusion take place at *T* = 100
  K, 200 K, 300 K and 500 K.

* For biological catalytic processes, a popular rule of thumb is
  that the rate doubles for every temperature increase of 10 K around
  room temperature.  What activation energy does this correspond to?

* Look at the relaxed configurations with the :command:`ase-gui`
  command::

    $ ase-gui -r 3,3,2 ontop.traj

  or::

    $ ase-gui -g 'd(4,-1),F[-1,2]' ontop.traj

  to plot the force in the *z*-direction on the gold atom as a
  function of the Au-Al distance. Note that -1 is the index of the last atom in the cell corresponding to the Au atom.  Try also *terminal-only-mode*::
 
    $ ase-gui -t -g 'd(4,-1),F[-1,2]' ontop.traj


Plot density differences
------------------------

It is sometimes useful to look at density changes when studying for
instance adsorption reactions. Copy the script
:download:`densitydiff.py` to your area.

Read it and try to understand what is does. Change the necessary lines
to look at one of your slabs with Au adsorbed. The script will write the
density difference to a :file:`.npy` file using NumPy's :func:`~numpy.save`
function (can be read with :func:`~numpy.load`.  Try this::
    
    from mayavi import mlab
    import numpy as np
    d = np.load('densitydiff.npy')
    d2 = np.tile(d, (2, 2, 1))  # repeat 2x2 times in x,y-plane
    mlab.contour3d(d2)
    mlab.show()
