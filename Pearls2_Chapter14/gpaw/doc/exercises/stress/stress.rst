.. _stress:

==========================================================
Plane wave mode and Stress tensor
==========================================================

The major advantage of running DFT calculations on a real space grid, is a very efficient paralleliation scheme when dealing with large systems. However, for small systems it is often faster to use a plane wave basis set instead. In this case all quantities are represented by their Fourier transforms on the periodic super cell and periodic boundary conditions are required. In grid mode, the key convergence parameter is the grid spacing `h`, whereas in planewave mode the corresponding parameter is `E_{cut}=G_{cut}^2/2` (in atomic units). `G_{cut}` determines the maximum size of reciprocal lattice vectors to be included in the plane wave expansion.  

Converging the plane wave cutoff
-------------------------------------

As a simple example we perform a calculation on bulk Si using the plane wave basis. The following script tests the convergence of the total energy with respect to the cutoff parameter `E_{cut}`.

.. literalinclude:: con_pw.py

Note that we use a rather rough k-point sampling and the default Fermi smearing of 0.1 eV. These parameters should of course be converged, but for now we will just keep them fixed in order to speed up the calculations. At `E_{cut}` = 400 eV, the total energy is seen to be converged to within 5 meV. Now edit the script such that the calculation is performed in grid mode with various values of the grid spacing (comment out the lines :samp:`for x in [100, ...` and :samp:`mode=PW(x)` and uncomment the lines :samp:`for h in [0.24, ...` and :samp:`h=x`). 

What grid spacing is needed in order to converge the total energy to within 5 meV? What is the timing of this calculation compared to the one at 400 eV in plane wave mode? Compare the number of grid points and plane waves used in these two calculations respectively (look in the output file 'convergence_400.txt' and the corresponding one for the grid calculation.)  You can see how long a calculation took by looking at the times recorded in the .txt file.

Optimizing the unit cell
------------------------

**Warning**: due to difficulties in optimizing cell and positions simultaneously
:class:`ase.constraints.UnitCellFilter` may produce
incorrect results. Always verify obtained structures by means of
performing separate cell
(see :class:`ase.constraints.StrainFilter`)
and positions optimizations (see :mod:`ase.optimize`).
Consider much more tighter fmax than the one used in this tutorial!

In the :ref:`aluminium_exercise` exercise the lattice constant of bulk Al was found by calculating the total energy at various lattice distances. A nice feature e of the plane wave mode is that it allows a simple implementation of the stress tensor, which can be used to optimize unit unit cells of periodic systems directly. The following script performs such an optimization for bulk Si.

.. literalinclude:: stress.py

The calculation uses 12 iterations to find the optimal lattice constant and the relaxation can be viewed with the command line tool ase-gui:

.. highlight:: bash

::

  $ ase-gui stress.txt

Since we know the experimental lattice constant, we could probably have calculated the PBE lattice constant faster by fitting a parabola to five points in the vicinity of the expermental lattice constant. However, for complicated unit cells with more than one lattice parameter, the stress tensor becomes a highly valuable tool.
