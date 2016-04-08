.. _mhtutorial:

================================================
Constrained minima hopping (global optimization)
================================================

This is an example of a search for a global optimum geometric configuration using the minima hopping algorithm, along with the Hookean class of constraints. This type of approach is useful in searching for the global optimum position of adsorbates on a surface while enforcing that the adsorbates' identity is preserved.

The below example looks at finding the optimum configuration of a :mol:`Cu_2` adsorbate on a fixed Pt (110) surface. Although this is not a physically relevant simulation --- these elements (Cu, Pt) were chosen only because they work with the EMT calculator -- one can imagine replacing the :mol:`Cu_2` adsorbate with CO, for example, to find its optimum binding configuration under the constraint that the CO does not dissociate into separate C and O adsorbates.

This also uses the Hookean constraint in two different ways. In the first, it constrains the Cu atoms to feel a restorative force if their interatomic distance exceeds 2.6 Angstroms; this preserves the dimer character of the :mol:`Cu_2`, and if they are near each other they feel no constraint. The second constrains one of the Cu atoms to feel a downward force if its position exceeds a z coordinate of 15 Angstroms. Since the Cu atoms are tied together, we don't necessarily need to put such a force on both of the Cu atoms. This second constraint prevents the :mol:`Cu_2` adsorbate from flying off the surface, which would lead to it exploring a lot of irrelevant configurational space, such as up in the vacuum or on the bottom of the next periodic slab.

.. literalinclude:: Cu2_Pt110.py

This script will produce 10 molecular dynamics and 11 optimization files. It will also produce a file called 'minima.traj' which contains all of the accepted minima. You can look at the progress of the algorithm in the file hop.log in combination with the trajectory files.

Alternatively, there is a utility to allow you to visualize the progress of the algorithm. You can run this from within the same directory as your algorithm as:

.. literalinclude:: mhsummary.py

This will make a summary figure, which should look something like the one below. As the search is inherently random, yours will look different than this (and this will look different each time the documentation is rebuilt). In this figure, you will see on the :math:`E_\mathrm{pot}` axes the energy levels of the conformers found. The flat bars represent the energy at the end of each local optimization step. The checkmark indicates the local minimum was accepted; red arrows indicate it was rejected for the three possible reasons. The black path between steps is the potential energy during the molecular dynamics (MD) portion of the step; the dashed line is the local optimization on termination of the MD step. Note the y axis is broken to allow different energy scales between the local minima and the space explored in the MD simulations. The :math:`T` and :math:`E_\mathrm{diff}` plots show the values of the self-adjusting parameters as the algorithm progresses.

.. image:: summary.png

You can see examples of the implementation of this for real adsorbates as well as find suitable parameters for the Hookean constraints:

  | Andrew Peterson
  | `Global optimization of adsorbate–surface structures while preserving molecular identity`__
  | Top. Catal., Vol. **57**, 40 (2014)

__ http://dx.doi.org/10.1007/s11244-013-0161-8

