.. _bandstructures:

=========================================
Calculation of electronic band structures
=========================================

In this tutorial we calculate the electronic band structure of Si along
high symmetry directions in the Brillouin zone.

First, a standard ground state calculations is performed and the results
are saved to a *.gpw* file. As we are dealing with small bulk system, 
plane wave mode is the most appropriate here.

.. literalinclude:: bandstructure.py
    :lines: 7-22

Next, :mod:`ase.dft.kpoints` module is used for generating k-points along
the high symmetry directions in the Brillouin zone. The below figure shows 
the high symmetry points of few common crystal lattices. 

.. figure:: ../../_static/bz-all.png
   :width: 600 px

For the band structure calculation, density is fixed to the previously
calculated ground state density (``fixdensity=True``), and as we want to
calculate all k-points, symmetry is not used (``symmetry='off'``). The
unoccupied states can be sometimes converged faster with the conjugate gradient
eigensolver.

.. literalinclude:: bandstructure.py
    :lines: 24-41

Finally, the bandstructure can be plotted e.g. with matplotlib. 
The :func:`ase.dft.kpoints.get_bandpath` provides in addition to the 
actual k-points information which is useful for plotting.

.. literalinclude:: bandstructure.py
    :lines: 43-61

.. figure:: bandstructure.png

The full script: :download:`bandstructure.py`.
