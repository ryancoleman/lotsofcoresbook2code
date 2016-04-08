.. _structure_optimizations:

======================
Structure optimization
======================
.. module:: ase.optimize
   :synopsis: Structure Optimization

The optimization algorithms can be roughly divided into local
optimization algorithms which find a nearby local minimum and
global optimization algorithms that try to find the global
minimum (a much harder task).


Local optimization
==================

The local optimization algorithms available in ASE are:
``BFGS``, ``LBFGS``, ``BFGSLineSearch``, ``LBFGSLineSearch``,
``MDMin``, and ``FIRE``.

.. seealso::

    `Performance test
    <https://wiki.fysik.dtu.dk/gpaw/devel/ase_optimize/ase_optimize.html>`_
    for all ASE local optimizers.


``MDMin`` and ``FIRE`` both use Newtonian dynamics with added
friction, to converge to an energy minimum, whereas the others are of
the quasi-Newton type, where the forces of consecutive steps are used
to dynamically update a Hessian describing the curvature of the
potential energy landscape.  You can use the ``QuasiNewton`` synonym
for ``BFGSLineSearch`` because this algorithm is in many cases the optimal
of the quasi-Newton algorithms.

All of the local optimizer classes have the following structure::

  class Optimizer:
      def __init__(self, atoms, restart=None, logfile=None):
      def run(self, fmax=0.05, steps=100000000):
      def get_number_of_steps():

The convergence criterion is that the force on all individual atoms
should be less than *fmax*:

.. math:: \max_a |\vec{F_a}| < f_\text{max}


BFGS
----
.. module:: ase.optimize.qn
   :synopsis: Quasi-Newton

The ``BFGS`` object is one of the minimizers in the ASE package. The below
script uses ``BFGS`` to optimize the structure of a water molecule, starting
with the experimental geometry::

   from ase import Atoms
   from ase.optimize import BFGS
   from ase.calculators.emt import EMT
   import numpy as np
   d = 0.9575
   t = np.pi / 180 * 104.51
   water = Atoms('H2O',
                 positions=[(d, 0, 0),
                            (d * np.cos(t), d * np.sin(t), 0),
                            (0, 0, 0)],
                 calculator=EMT())
   dyn = BFGS(water)
   dyn.run(fmax=0.05)

which produces the following output. The columns are the solver name, step
number, clock time, potential energy (eV), and maximum force.::

   BFGS:   0  19:45:25        2.769633       8.6091
   BFGS:   1  19:45:25        2.154560       4.4644
   BFGS:   2  19:45:25        1.906812       1.3097
   BFGS:   3  19:45:25        1.880255       0.2056
   BFGS:   4  19:45:25        1.879488       0.0205


When doing structure optimization, it is useful to write the
trajectory to a file, so that the progress of the optimization run can
be followed during or after the run::

  dyn = BFGS(water, trajectory='H2O.traj')
  dyn.run(fmax=0.05)

Use the command ``ase-gui H2O.traj`` to see what is going on (more here:
:mod:`ase.gui`).  The trajectory file can also be accessed using the
module :mod:`ase.io.trajectory`.

The ``attach`` method takes an optional argument ``interval=n`` that can
be used to tell the structure optimizer object to write the
configuration to the trajectory file only every ``n`` steps.

During a structure optimization, the BFGS and LBFGS optimizers use two
quantities to decide where to move the atoms on each step:

* the forces on each atom, as returned by the associated
  :class:`~ase.calculators.calculator.Calculator` object
* the Hessian matrix, i.e. the matrix of second derivatives
  :math:`\frac{\partial^2 E}{\partial x_i \partial x_j}` of the
  total energy with respect to nuclear coordinates.

If the atoms are close to the minimum, such that the potential energy
surface is locally quadratic, the Hessian and forces accurately
determine the required step to reach the optimal structure.  The
Hessian is very expensive to calculate *a priori*, so instead the
algorithm estimates it by means of an initial guess which is adjusted
along the way depending on the information obtained on each step of
the structure optimization.

It is frequently practical to restart or continue a structure
optimization with a geometry obtained from a previous relaxation.
Aside from the geometry, the Hessian of the previous run can and
should be retained for the second run.  Use the ``restart`` keyword to
specify a file in which to save the Hessian::

  dyn = BFGS(atoms=system, trajectory='qn.traj', restart='qn.pckl')

This will create an optimizer which saves the Hessian to
:file:`qn.pckl` (using the Python :mod:`pickle` module) on each
step.  If the file already exists, the Hessian will also be
*initialized* from that file.

The trajectory file can also be used to restart a structure
optimization, since it contains the history of all forces and
positions, and thus whichever information about the Hessian was
assembled so far::

  dyn = BFGS(atoms=system, trajectory='qn.traj')
  dyn.replay_trajectory('history.traj')

This will read through each iteration stored in :file:`history.traj`,
performing adjustments to the Hessian as appropriate.  Note that these
steps will not be written to :file:`qn.traj`.  If restarting with more than
one previous trajectory file, use :command:`ase-gui` to concatenate them
into a single trajectory file first::

  $ ase-gui part1.traj part2.traj -o history.traj

The file :file:`history.traj` will then contain all necessary
information.

When switching between different types of optimizers, e.g. between
``BFGS`` and ``LBFGS``, the pickle-files specified by the
``restart`` keyword are not compatible, but the Hessian can still be
retained by replaying the trajectory as above.


LBFGS
-----
.. module:: ase.optimize.lbfgs

LBFGS is the limited memory version of the BFGS algorithm, where
the inverse of Hessian matrix is updated instead of the Hessian
itself. Two ways exist for determining the atomic
step: Standard ``LBFGS`` and ``LBFGSLineSearch``. For the
first one, both the directions and lengths of the atomic steps
are determined by the approximated Hessian matrix. While for the
latter one, the approximated Hessian matrix is only used to find
out the directions of the line searches and atomic steps, the
step lengths are determined by the forces.

To start a structure optimization with LBFGS algorithm is similar to
BFGS. A typical optimization should look like::

  dyn = LBFGS(atoms=system, trajectory='lbfgs.traj', restart='lbfgs.pckl')

where the trajectory and the restart save the trajectory of the
optimization and the vectors needed to generate the Hessian Matrix.


FIRE
----
.. module:: ase.optimize.fire

Read about this algorithm here:

  | Erik Bitzek, Pekka Koskinen, Franz Gähler, Michael Moseler, and Peter Gumbsch
  | `Structural Relaxation Made Simple`__
  | Physical Review Letters, Vol. **97**, 170201 (2006)

__ http://dx.doi.org/10.1103/PhysRevLett.97.170201


MDMin
-----
.. module:: ase.optimize.mdmin

The MDmin algorithm is a modification of the usual velocity-Verlet
molecular dynamics algorithm.  Newtons second law is solved
numerically, but after each time step the dot product between the
forces and the momenta is checked.  If it is zero, the system has just
passed through a (local) minimum in the potential energy, the kinetic
energy is large and about to decrease again.  At this point, the
momentum is set to zero.  Unlike a "real" molecular dynamics, the
masses of the atoms are not used, instead all masses are set to one.

The MDmin algorithm exists in two flavors, one where each atom is
tested and stopped individually, and one where all coordinates are
treated as one long vector, and all momenta are set to zero if the
dot product between the momentum vector and force vector (both of
length 3N) is zero.  This module implements the latter version.

Although the algorithm is primitive, it performs very well because it
takes advantage of the physics of the problem.  Once the system is so
near the minimum that the potential energy surface is approximately
quadratic it becomes advantageous to switch to a minimization method
with quadratic convergence, such as *Conjugate Gradient* or *Quasi
Newton*.


SciPy optimizers
----------------
.. module:: ase.optimize.sciopt

SciPy provides a number of optimizers. An interface module for a couple of
these have been written for ASE. Most notable are the optimizers SciPyFminBFGS
and SciPyFminCG. These are called with the regular syntax and can be imported
as::

  from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG

.. autoclass:: ase.optimize.sciopt.SciPyFminBFGS
.. autoclass:: ase.optimize.sciopt.SciPyFminCG

.. seealso::

  :epydoc:`optimize.sciopt.SciPyFminBFGS`,
  :epydoc:`optimize.sciopt.SciPyFminCG`


BFGSLineSearch
--------------
.. module:: ase.optimize.bfgslinesearch

BFGSLineSearch is the BFGS algorithm with a line search mechanism
that enforces the step taken fulfills the Wolfe conditions, so that
the energy and absolute value of the force decrease monotonically. Like
the LBFGS algorithm the inverse of the Hessian Matrix is updated.

The usage of BFGSLineSearch algorithm is similar to other BFGS type
algorithms. A typical optimization should look like::

  from ase.optimize.bfgslinesearch import BFGSLineSearch

  dyn = BFGSLineSearch(atoms=system, trajectory='bfgs_ls.traj', restart='bfgs_ls.pckl')

where the trajectory and the restart save the trajectory of the
optimization and the information needed to generate the Hessian Matrix.

.. note::

   In many of the examples, tests, exercises and tutorials,
   ``QuasiNewton`` is used -- it is a synonym for ``BFGSLineSearch``.

The BFGSLineSearch algorithm is not compatible with nudged elastic band
calculations.


Global optimization
===================

There are currently two global optimisation algorithms available.


Basin hopping
-------------
.. module:: ase.optimize.basin

The global optimization algorithm can be used quite similar as a
local optimization algorithm::

  from ase import *
  from ase.optimize.basin import BasinHopping

  bh = BasinHopping(atoms=system,         # the system to optimize
                    temperature=100 * kB, # 'temperature' to overcome barriers
                    dr=0.5,               # maximal stepwidth
                    optimizer=LBFGS,      # optimizer to find local minima
                    fmax=0.1,             # maximal force for the optimizer
                    )

Read more about this algorithm here:

  | David J. Wales and Jonathan P. K. Doye
  | `Global Optimization by Basin-Hopping and the Lowest Energy Structures of Lennard-Jones Clusters Containing up to 110 Atoms`__
  | J. Phys. Chem. A, Vol. **101**, 5111-5116 (1997)

__ http://pubs.acs.org/doi/abs/10.1021/jp970984n

and here:

  | David J. Wales and Harold A. Scheraga
  | `Global Optimization of Clusters, Crystals, and Biomolecules`__
  | Science, Vol. **285**, 1368 (1999)

__ http://www.sciencemag.org/cgi/content/abstract/sci;285/5432/1368

Minima hopping
--------------

The minima hopping algorithm was developed and described by Goedecker:

  | Stefan Goedecker
  | `Minima hopping: An efficient search method for the global minimum of the potential energy surface of complex molecular systems`__
  | J. Chem. Phys., Vol. **120**, 9911 (2004)

__ http://dx.doi.org/10.1063/1.1724816

This algorithm utilizes a series of alternating steps of NVE molecular dynamics and local optimizations, and has two parameters that the code dynamically adjusts in response to the progress of the search. The first parameter is the initial temperature of the NVE simulation. Whenever a step finds a new minimum this temperature is decreased; if the step finds a previously found minimum the temperature is increased. The second dynamically adjusted parameter is :math:`E_\mathrm{diff}`, which is an energy threshold for accepting a newly found minimum. If the new minimum is no more than :math:`E_\mathrm{diff}` eV higher than the previous minimum, it is acccepted and :math:`E_\mathrm{diff}` is decreased; if it is more than :math:`E_\mathrm{diff}` eV higher it is rejected and :math:`E_\mathrm{diff}` is increased. The method is used as::

   from ase.optimize.minimahopping import MinimaHopping
   opt = MinimaHopping(atoms=system)
   opt(totalsteps=10)

This will run the algorithm until 10 steps are taken; alternatively, if totalsteps is not specified the algorithm will run indefinitely (or until stopped by a batch system). A number of optional arguments can be fed when initializing the algorithm as keyword pairs. The keywords and default values are:


 | ``T0``: 1000.,  # K, initial MD 'temperature'
 | ``beta1``: 1.1,  # temperature adjustment parameter
 | ``beta2``: 1.1,  # temperature adjustment parameter
 | ``beta3``: 1. / 1.1,  # temperature adjustment parameter
 | ``Ediff0``: 0.5,  # eV, initial energy acceptance threshold
 | ``alpha1`` : 0.98,  # energy threshold adjustment parameter
 | ``alpha2`` : 1. / 0.98,  # energy threshold adjustment parameter
 | ``mdmin`` : 2,  # criteria to stop MD simulation (no. of minima)
 | ``logfile``: 'hop.log',  # text log
 | ``minima_threshold`` : 0.5,  # A, threshold for identical configs
 | ``timestep`` : 1.0,  # fs, timestep for MD simulations
 | ``optimizer`` : QuasiNewton,  # local optimizer to use
 | ``minima_traj`` : 'minima.traj',  # storage file for minima list

Specific definitions of the ``alpha``, ``beta``, and ``mdmin`` parameters can be found in the publication by Goedecker. ``minima_threshold`` is used to determine if two atomic configurations are identical; if any atom has moved by more than this amount it is considered a new configuration. Note that the code tries to do this in an intelligent manner: atoms are considered to be indistinguishable, and translations are allowed in the directions of the periodic boundary conditions. Therefore, if a CO is adsorbed in an ontop site on a (211) surface it will be considered identical no matter which ontop site it occupies.

The trajectory file ``minima_traj`` will be populated with the accepted minima as they are found. A log of the progress is kept in ``logfile``.

The code is written such that a stopped simulation (e.g., killed by the batching system when the maximum wall time was exceeded) can usually be restarted without too much effort by the user. In most cases, the script can be resubmitted without any modification -- if the ``logfile`` and ``minima_traj`` are found, the script will attempt to use these to resume. (Note that you may need to clean up files left in the directory by the calculator, however, such as the .nc file produced by Jacapo.)

Note that these searches can be quite slow, so it can pay to have multiple searches running at a time. Multiple searches can run in parallel and share one list of minima. (Run each script from a separate directory but specify the location to the same absolute location for ``minima_traj``). Each search will use the global information of the list of minima, but will keep its own local information of the initial temperature and :math:`E_\mathrm{diff}`.

For an example of use, see the :ref:`mhtutorial` tutorial.
