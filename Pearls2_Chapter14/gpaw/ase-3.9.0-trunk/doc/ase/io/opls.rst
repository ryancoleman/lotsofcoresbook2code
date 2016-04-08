==========================================
Setting up an OPLS force field calculation
==========================================

.. module:: ase.io.opls
   :synopsis: OPLS force field

In order to facilitate the definition of structures for the use
of OPLS force fields, there are some helper classes.


Modified xyz
============

Suppose, we define the ethanal molecule as an modified xyz file 
(``172_mod.xyz``):

.. literalinclude:: 172_mod.xyz

Then we can read and view the structure using:

.. literalinclude:: view_172_mod.py


Defining the force field
========================

The definitions of the force field can be stored in an Amber like style 
(``172_defs.par``):

.. literalinclude:: 172_defs.par

We can write LAMMPS input using the information above:

.. literalinclude:: write_lammps.py

which writes the LAMMPS input files ``lmp_atoms`` defining atoms,  bonds,
etc., and ``lmp_opls`` defining the corresponding OPLS force field. A
rudimentary ``lmp_in`` is also written.
