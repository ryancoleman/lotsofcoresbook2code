.. _fcc_alloys_tutorial:

===============================
GA Search for stable FCC alloys
===============================

In this tutorial we will emulate an older paper [Jóhannesson]_ and determine
the most stable FCC alloy using the genetic algorithm. Since the purpose is
only the tutorial we will limit the phase space to the elements supported by
the `EMT potential`_. The search is also equivalent to the recent search for
mixed metal ammines with superior properties for ammonia storage described
here:

.. _`EMT potential`: https://wiki.fysik.dtu.dk/ase/ase/calculators/emt.html#module-ase.calculators.emt

   | P. B. Jensen, S. Lysgaard, U. J. Quaade and T. Vegge
   | `Designing Mixed Metal Halide Ammines for Ammonia Storage Using Density Functional Theory and Genetic Algorithms`__
   | Physical Chemistry Chemical Physics, Vol **16**, No. 36, pp. 19732-19740, (2014)
   
   __ http://dx.doi.org/10.1039/C4CP03133D

.. contents::
   
   
Basic outline of the search
===========================

1. Choose the phase space of your problem. Is the number of possible
   individuals large enough to prevent a full screening and is the fitness
   function too discontinuous for a traditional optimization by derivation? If
   so continue.

2. Choose model structures and calculate references in those structures. Put
   the results somewhere accesible for a script initiated by the genetic
   algorithm.

3. Choose suitable parameters like population size (general rule of thumb for
   the population size: `log_2(N)` < pop size < `2log_2(N)`, where `N` is the
   size of the phase space), convergence criteria etc.

4. Create the initial population.

5. Choose procreation operators, i.e. how should offspring be produced. New
   operators can easily be created by modifying the existing operators.

6. Run the algorithm.

Here we would like to predict the most stable fcc alloys. In this tutorial we
only have the `EMT potential`__ available thus we are limited to the
supported metal elements: Al, Ni, Cu, Pd, Ag, Pt and Au. We limit ourselves
to at most 4 different metals in one structure, thereby having only `7^4 =
2401` candidates in the phase space, symmetry would make this number even
lower but the number is fitting for this tutorial.

__ `EMT potential`_

For a real application of the algorithm it is necessary to use a more sophisticated calculator, in that case each individual calculation is performed on a cluster by submitting to a queuing system. How this is achieved in the algorithm is covered in :ref:`genetic_algorithm_optimization_tutorial`.

.. defined for an alloy :mol:`ABC_2`: A + B + 2C -> :mol:`ABC_2` as: `\Delta H_f = E_{ABC2} - E_A - E_B - 2E_C`


.. _references:

Setting up reference database
=============================

Now we need to set up a database in which
reference calculations can be stored. This can either
be in a central database server where keywords distinguish
between different references or dedicated separate
databases for each different type of reference calculations.

In the following script, :download:`ga_fcc_references.py`, we put the references in the database file *refs.db*. Our model structure is fcc which is loaded with :func:`ase.lattice.cubic.FaceCenteredCubic`. We perform a volume relaxation to find the optimal lattice constant and lowest energy, which we save in the database as key-value pairs for quick retrieval.

.. literalinclude:: ga_fcc_references.py
                                        
                                        
Initial population
==================

We choose a population size of 10 individuals and create the initial population by randomly selecting four elements for each starting individual.

.. literalinclude:: ga_fcc_alloys_start.py
                                        
Note how we add the population size and metals as extra key-value pairs when we create the database *fcc_alloys.db*. We can then retrieve these parameters later when running the main script to avoid having to input the same parameters twice.

We can study our initial population by doing (on the command-line)::
  
    $ ase-db fcc_alloys.db -c +atoms_string
        
the term ``atoms_string`` determines the order in which the elements are put into the model structure. So it is possible to fully describe an individual by just providing the ``atoms_string``.
        
        
.. _`main script`:

Run the algorithm
=================

.. literalinclude:: ga_fcc_alloys_main.py
                                        
In this script we run a generational GA as opposed to the pool GA outlined in :ref:`genetic_algorithm_optimization_tutorial`. This is achieved by having two for-loops; the innermost loop runs the number of times specified by the population size it corresponds to one generation. The outermost loop runs as many generations as specified in ``num_gens``. The function :func:`pop.update()` is called after the innermost loop has finished thereby only adding individuals to the population after a whole generation is calculated.

After each generation is finished the population is printed to the screen so we can follow the evolution. The calculated individuals are continuously added to ``fcc_alloys.db``, we can evaluate them directly by doing from the command line (in another shell instance if the GA is still running)::

    $ ase-db fcc_alloys.db -c +atoms_string,raw_score,generation,hof -s raw_score

*Note:* When reading the database using ase-db, it might be necessary to increase the number of shown entries, e.g. ``ase-db fcc-alloys.db --limit N``, where ``N`` is the number of entries to show (as default the first 500 entries are shown, ``--limit 0`` will show all. For further info use the help: ``ase-db –help``, or consult the `ase-db manual`_).

.. _`ase-db manual`: https://wiki.fysik.dtu.dk/ase/ase/db/db.html#module-ase.db

To prevent clutter we import the relax function from the following script:

.. _`relaxation script`:
                                        
.. literalinclude:: ga_fcc_alloys_relax.py
                                        
The relaxation script is naturally similar to the script we used to calculate the references_.

*Note* that the global optimum is :mol:`PtNi_3` with a -0.12 eV heat of formation, whereas the second worst alloy is :mol:`AlNi_3` heat of formation 0.26 eV. This result is in complete contrast to the conclusion obtained in [Jóhannesson]_, where :mol:`AlNi_3` is the most stable alloy within the phase space chosen here. Obviously there is a limit to the predictive power of EMT!
                                        
Extending the algorithm
=======================

There are different ways one can extend the algorithm and make it more complex and sophisticated, all employed in [Jensen]_:


Extra mutation operators
------------------------

Instead of only using random operations we can include some that mutates elements to other elements nearby in the periodic table::

  from ase.ga.element_mutations import RandomElementMutation
  from ase.ga.element_mutations import MoveDownMutation
  from ase.ga.element_mutations import MoveUpMutation
  from ase.ga.element_mutations import MoveLeftMutation
  from ase.ga.element_mutations import MoveRightMutation
  from ase.ga.element_crossovers import OnePointElementCrossover
  
  ...
  
  oclist = ([4,1,1,1,1,8], [RandomElementMutation([metals]),
                            MoveDownMutation([metals]),
                            MoveUpMutation([metals]),
                            MoveLeftMutation([metals]),
                            MoveRightMutation([metals]),
                            OnePointElementCrossover([metals])])
  mut_selector = MutationSelector(*oclist)

These operators takes advantage of the fact that chemically like elements (close in the periodic table) exhibit similar properties and the substitution of one to a chemically similar elements could refine the properties of an alloy in the population. A natural extension of these operators would be to use a different ordering of the elements than the periodic table; e.g. Pettifor chemical scale, electronegativity, etc.

Note how we have set the probabilities for selecting operators differently. The probability for ``RandomElementMutation`` is equal to the sum of the *move* mutations. Similarly the probability of ``OnePointElementCrossover`` is equal to the sum of all the mutation operators. This is to prevent the search from being purely local.


Prevent identical calculations from being performed
---------------------------------------------------

In the current `main script`_ there is no check to determine whether an identical calculation has been performed, this is easy to check in this regime where model structures are used and we can just use the ``atoms_string``. We insert the following in the inner loop::

  for i in range(population_size):
      dup = True
      while dup:
          a1, a2 = pop.get_two_candidates(with_history=False)
          op = operation_selector.get_operator()
          a3, desc = op.get_new_individual([a1, a2])

          dup = db.is_duplicate(atoms_string=''.join(a3.get_chemical_symbols()))

Since the fcc model structure is completely symmetric we could compare sorted versions of the ``atoms_string``, thereby ruling out individuals containing the same elements in different order.


Reuse of calculations between algorithm runs
--------------------------------------------

Since genetic algorithms are inherently random in nature one can never be sure to obtain the global minimum with only one algorithm run, it is customary to perform more runs and check that the results agree. In this case it is vital to be able to reuse identical calculations between runs.

We do the following from the command line to create a new database file containing only the relaxed structures::

    $ ase-db fcc_alloys.db relaxed=1 -i all_relaxed.db
        
We subsequently add this to the `relaxation script`_::
  
  def relax(input_atoms, ref_db):
      atoms_string = input_atoms.get_chemical_symbols()
      relaxed_db = connect('all_relaxed.db')
      save_relax = True
      try:
          dct = relaxed_db.get(atoms_string=''.join(atoms_string))
      except KeyError:
          # Open connection to the database with reference data
          db = connect(ref_db)
                  
      # Omitting lines up to the point where hof has been calculated
          ...
          
      else:
          hof = dct.hof
          latticeconstant = dct.latticeconstant
          save_relax = False
      # Place the calculated parameters in the info dictionary of the
      # input_atoms object
          
      ...
          
      # Put this at the very end
      if save_relax:
          relaxed_db.write(input_atoms,relaxed=1,
                           key_value_pairs=input_atoms.info['key_value_pairs'])
  
Before the actual calculation is performed ``all_relaxed.db`` is checked to see if it has been calculated before; if so we just collect the heat of formation, but if not we do the calculation and save it directly to ``all_relaxed.db``.
*Note:* this addition assumes that `Prevent identical calculations from being performed`_.

.. [Jóhannesson] G. Jóhannesson, T. Bligaard, A. Ruban, H. Skriver, K. Jacobsen and J. Nørskov.
   Combined Electronic Structure and Evolutionary Search Approach to Materials Design,
   Phys. Rev. Lett., Vol **88**, No. 25, pp. 1-5 (2002)
.. [Jensen] P. B. Jensen, S. Lysgaard, U. J. Quaade and T. Vegge.
   Designing Mixed Metal Halide Ammines for Ammonia Storage Using Density Functional Theory and Genetic Algorithms
   Phys. Chem. Chem. Phys., Vol **16**, No. 36, pp. 19732-19740, (2014)
