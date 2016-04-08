.. _plot_wave_functions:

=======================
Plotting wave functions
=======================

-----------------------------
Creating a wave function file
-----------------------------

The following script will do a calculation for a CO
molecule and save the wave functions in a file (``CO.gpw``).

.. literalinclude:: CO.py

---------------------------------
Creating wave function cube files
---------------------------------

You can get seperate cube files (the format used by Gaussian) for each wavefunction with the script:

.. literalinclude:: CO2cube.py

The script produced the files CO_0.cube .. CO_5.cube, which might be viewed using for example `jmol <http://jmol.sourceforge.net/>`_ or `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_. 

Plotting wave functions with jmol
---------------------------------

See http://jmol.sourceforge.net/docs/surface/

You can visualize an isosurface from a cube file by creating the
myscript.spt script::
  
   #isosurface name cutoff "cubefile"
   isosurface pos 0.05 "CO_1.cube"

and executing it with jmol::

   jmol -s myscript.spt

You can also save the image directly on the command line::

   jmol -ions myscript.spt -w JPEG:myfile.jpg
  
Plotting wave functions with VMD
--------------------------------

To view the wavefunctions, start VMD with the command line::

    vmd CO_*.cube

You will get two windows, one with a very primitive representation of
the molecule, and one called 'VMD Main'.  In the VMD Main window,
select Graphics/Representation.  In the top panel there will be a
single "representation" of the molecule called "Line".  Select it, and
change the type to "CPK".  This will show the molecule with spheres for
atoms and rods for bonds.  You probably want to reduce the size of the
spheres and increase the resolution of the spheres.  Now you see the
atoms!

To see the wavefunctions, click on the "Create Rep" button, and select
the new representation.  Select type "Isosurface".  Near the bottom,
find the menu labeled "Draw" and select the value "Solid Surface".
Now you can see an iso-surface of the wavefunction, you select the
value in the field "Isovalue".  The default 0 is rarely useful.  You
select the different Kohn-Sham states (stored in the different cube
files) with the pull-down menu labeled "Vol".

**IMPORTANT:** This works best for molecules.  In solids, the
wavefunctions will be complex, VMD does not handle this well.

Creating cube to plt files (gOpenMol)
-------------------------------------

**Warning** this section is obsolete!

The cube files can be transformed to plt format using the program g94cub2pl from the gOpenMol utilities.

----------------------------------------------
Creating wave function plt files with gOpenMol
----------------------------------------------

**Warning** this section is obsolete!

One can write out the wave functions in the very compact (binary) `gOpenMol <http://www.csc.fi/gopenmol/>`_ plt format directly:

.. literalinclude:: CO2plt.py

