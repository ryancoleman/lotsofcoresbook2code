.. _aluminium_exercise:

=============
Bulk aluminum
=============

We will look at bulk fcc aluminum and make a single energy calculation at the
experimental lattice constant `a_0` = 4.05 Å. For the first example, we
choose a plane-wave cutoff energy of 300 eV and 4 x 4 x 4 **k**-points.  Copy
this :download:`Al_fcc.py` to a place in your file area:

.. literalinclude:: Al_fcc.py

.. highlight:: bash

Read the script and try to get an idea of what it will do. Run the
script by typing::

  $ python Al_fcc.py

The program will pop up a window showing the bulk structure.  Verify
that the structure indeed is fcc. Try to identify the closepacked (111)
planes.  In :program:`ase-gui` this is done by choosing :menuselection:`View
--> Repeat`.

Notice that the program has generated two output files::

  Al-fcc.gpw
  Al-fcc.txt

* A tar-file (conventional suffix :file:`.gpw`) containing binary data
  such as eigenvalues, electron density and wave functions (see
  :ref:`restart_files`).

* An ASCII formatted log file (conventional suffix :file:`.txt`) that
  monitors the progress of the calculation.

Try to take a look at the file :file:`Al-fcc.txt`.  Find the number of
plane-waves used.  You can conveniently monitor some variables by using the
:command:`grep` utility.  By typing::

  $ grep iter Al-fcc.txt

you see the progress of the iteration cycles including convergence of
wave functions, density and total energy. If the ``txt`` keyword is omitted 
the log output will be printed directly in the terminal.

.. highlight:: python

The binary file contains all information about the calculation. Try
typing the following from the Python interpreter:

>>> from gpaw import GPAW
>>> calc = GPAW('Al-fcc.gpw', txt=None)
>>> bulk = calc.get_atoms()
>>> print bulk.get_potential_energy()
-4.12234332252
>>> density = calc.get_pseudo_density()
>>> density.shape
(9, 9, 9)
>>> density.max()
0.42718359271458561
>>> from mayavi import mlab
>>> mlab.contour3d(density)
<mayavi.modules.iso_surface.IsoSurface object at 0x7f1194491110>
>>> mlab.show()


Equilibrium lattice properties
==============================

We now proceed to calculate some equilibrium lattice properties of
bulk Aluminum.

* First map out the cohesive curve `E(a)` for Al(fcc), i.e.  the
  total energy as function of lattice constant `a`, around the
  experimental equilibrium value of `a_0` = 4.05 Å.  Get four or more
  energy points, so that you can make a fit.

  .. hint::

     Modify :download:`Al_fcc.py` by adding a
     for-loop like this::

         for a in [3.9, 4.0, 4.1, 4.2]:
             name = 'bulk-fcc-%.1f' % a

     and then indent the rest of the code by
     four spaces (Python uses indentation to group statements together
     - thus the for-loop will end at the first unindented line).
     Remove the ``view(bulk)`` line.

* Fit the data you have obtained to get `a_0` and the energy curve
  minimum `E_0=E(a_0)`.  From your fit, calculate the bulk
  modulus

  .. math:: B = V\frac{d^2 E}{dV^2} = \frac{M}{9a_0}\frac{d^2 E}{da^2},

  where *M* is the number of atoms per cubic unit cell:
  `V=a^3/M` (`M=4` for fcc).  Make the fit using your favorite math
  package (NumPy_) or use :program:`ase-gui` like this::

    $ ase-gui bulk-*.txt

  Then choose :menuselection:`Tools --> Bulk Modulus`.

  Another alternative is to use the :ase:`Equation of state module
  <ase/utils.html#equation-of-state>` (see this :ase:`tutorial
  <tutorials/eos/eos.html>`).

* Compare your results to the experimental values `a_0` = 4.05 Å and `B`
  = 76 GPa.  Mind the units when you calculate the bulk modulus (read
  about ASE-units :ase:`here <ase/units.html>`).
  What are the possible error sources?

  .. note::

      The LDA reference values are: `a_0` = 3.98 Å and `B` = 84.0 GPa -
      see S. Kurth *et al.*, Int. J. Quant. Chem. **75** 889-909
      (1999).

.. _NumPy: http://docs.scipy.org/doc/numpy/reference/generated/
    numpy.polynomial.polynomial.Polynomial.html
    #numpy-polynomial-polynomial-polynomial

    
Convergence in number of **k**-points
-------------------------------------

Now we will investigate the necessary **k**-point sampling for bulk
fcc Aluminum; this is a standard first step in all DFT calculations.

* Repeat the calculation above for the equilibrium lattice constant
  for more dense Brillouin zone samplings (try ``k=6,8,10,...``).

  .. hint::

      You may want to speed up these calculations by running them in
      parallel.

* Estimate the necessary number of **k**-points for achieving an
  accurate value for the lattice constant.

* Do you expect that this **k**-point test is universal for all other
  Aluminum structures than fcc?  What about other chemical elements ?


Equilibrium lattice properties for bcc
======================================

* Set up a similar calculation for bcc, in the minimal unit cell.
  
* Make a qualified starting guess on `a_\text{bcc}` from the lattice
  constant for fcc, that you have determined above. One can either
  assume that the primitive unit cell volumes of the fcc and bcc
  structure are the same or that the nearest neighbor distances are
  the same. Find a guess for `a_\text{bcc}` for both
  assumptions. Later, you can comment on which assumption gives the
  guess closer to the right lattice constant.

* Check that your structure is right by repeating the unit cell.

* Map out the cohesive curve `E(a)` for Al(bcc) and determine
  `a_\text{bcc}`, using a few points.  Is it a good idea to use the
  same **k**-point setup parameters as for the fcc calculations?
  Calculate the bulk modulus, as it was done for fcc, and compare the
  result to the fcc bulk modulus. What would you expect?

* Using the lattice constants determined above for fcc and bcc,
  calculate the fcc/bcc total energies.  The total energies that GPAW
  calculates are relative to isolated atoms (more details here:
  :ref:`zero_energy`).  This exercise is sensitive to the number of
  **k**-points, make sure that your **k**-point sampling is dense
  enough.  Also make sure your energies are converged with respect to
  the plane-wave cutoff energy.
