.. _wannier:
    
=================
Wannier Functions
=================

In order to get a feel for chemical bonds in molecules and solids, we
can transform the Kohn-Sham orbitals into a set of maximally localized
Wannier functions.  We have cheated a little bit and prepared a file
for bulk silicon and a Benzene molecule so that you only have to
concentrate on the wannier analysis of the molecules.

Start by running :svn:`~doc/exercises/wannier/si.py` and
make sure you agree with the way the diamond structure is set up. The
resulting :file:`.gpw` file is used as input to
:svn:`~doc/exercises/wannier/wannier-si.py` which
transforms the Kohn-Sham orbitals to maximally localized wannier
functions and plot the atoms along with the centers of the wannier
functions.  Note that the wannier centers are treated as "X" atoms
which are plotted as small red spheres.  How many covalent bonds do
you expect in a unit cell with 8 tetravalent Silicon atoms?

The script :svn:`~doc/exercises/wannier/benzene.py`
produces a :file:`.gpw` that can be used as input to create wannier
functions. Convince yourself that the chosen number of bands matches
the number of occupied orbitals in the molecule.  How many covalent
bonds do you expect in Benzene?  Look at
:svn:`~doc/exercises/wannier/wannier-benzene.py` and figure
out what it does. Run it and look at the graphical representation.
Note in particular the alternating single/double bonds between the
carbon atoms.  What happens if also you include one or two unoccupied
bands?  The script also produces two :file:`.cube` files. One contains
the wavefunction of the Highest Occupied Molecular Orbital (HOMO) and
the other contains a wannier function centered between a Carbon and a
Hydrogen atom. Study these with :program:`VMD` and determine which
type of orbitals they represent (:math:`\sigma` or :math:`\pi`).

Now repeat the wannier function analysis on the following molecules

* H2O : use your own files from the vibrational exercise, but make
  sure the number of bands is equal to the number of occupied orbitals.

* CO : use your own :file:`.gpw` file from the wavefunction
  exercise. Is it a single, double or triple bond?

or study your own favorite molecule.

.. hint::
  
  To be able to see the Wannier centers, it might be necessary to
  decrease the atomic radii, so the spheres don't overlap.
  In :program:`ase-gui` this can be done by choosing 
  :menuselection:`View --> Settings`, and
  then decrease the scaling factor of the covalent radii.
