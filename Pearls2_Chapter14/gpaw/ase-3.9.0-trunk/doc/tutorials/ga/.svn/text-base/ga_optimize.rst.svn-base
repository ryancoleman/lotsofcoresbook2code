.. _genetic_algorithm_optimization_tutorial:

=====================================
Optimization with a Genetic Algorithm
=====================================

A genetic algorithm (GA) has been implemented for global structure
optimization within ase. The optimizer consists of its own module
:mod:`ase.ga` which includes all classes needed for the optimizer.

The method was first described in the supplemental material of

   | L. B. Vilhelmsen and B. Hammer
   | `Systematic Study of Au6 to Au12 Gold Clusters on MgO(100) F Centers Using Density-Functional Theory`__
   | Physical Review Letters, Vol. 108 (Mar 2012), 126101

   __ http://dx.doi.org/10.1103/physrevlett.108.126101

and a full account of the method is given in

   | L. B. Vilhelmsen and B. Hammer
   | `A genetic algorithm for first principles global optimization of supported nano structures`__
   | Journal of Chemical Physics, Vol 141, 044711 (2014)

   __ http://dx.doi.org/10.1063/1.4886337

Any questions about how to use the GA can be asked at the mailing
list.


A Brief Overview of the Implementation
======================================

The GA relies on the ase.db module for tracking which structures have
been found. Before the GA optimization starts the user therefore needs
to prepare this database and appropriate folders. This is done trough
an initialization script as the one described in the next section. In
this initialization the starting population is generated and
added to the database.

After initialization the main script is run. This script defines
objects responsible for the different parts of the GA and then creates
and locally relaxes new candidates. It is up to the user to define
when the main script should terminate. An example of a main script is
given in the next section.  Notice that because of the persistent data
storage the main script can be executed multiple times to generate new
candidates.

The GA implementation generally follows a responsibility driven
approach. This means that each part of the GA is isolated into
individual classes making it possible to put together an optimizer
satisfying the needs of a specific optimization problem.

This tutorial will use the following parts of the GA:

* A population responsible for proposing new candidates to pair
  together.
* A paring operator which combines two candidates.
* A set of mutations.
* A comparator which determines if two structures are different.
* A starting population generator.

Each of the above components are described in the supplemental
material of the first reference given above and will not be discussed
here. The example will instead focus on the technical aspect of
executing the GA.

A Basic Example
===============
The user needs to specify the following three properties about the
structure that needs to be optimized.

* A list of atomic numbers for the structure to be optimized

* A super cell in which to do the optimization. If the structure to
  optimize resides on a surface or in a support this supercell
  contains the atoms which should not be considered explicitly by the
  GA.

* A box defining the volume of the super cell in which to randomly
  distribute the starting population.

As an example we will find the structure of a
:mol:`Ag_2Au_2` cluster on a Au(111) surface using the
EMT optimizer.

The script doing all the initialisations should be run in the folder
in which the GA optimisation is to take place. The script looks as follows:

.. literalinclude:: basic_example_create_database.py

Having initialized the GA optimization we now need to actually run the
GA. The main script running the GA consists of first an initialization
part, and then a loop proposing new structures and locally optimizing
them. The main script can look as follows:

.. literalinclude:: basic_example_main_run.py

The above script proposes and locally relaxes 20 new candidates. To
speed up the execution of this sample the local relaxations are
limited to 100 steps. This restriction should not be set in a real
application. *Note* it is important to set the the ``raw_score``, as
it is what is being optimized (maximized). It is really an input in the
``atoms.info['key_value_pairs']`` dictionary.

The GA progress can be monitored by running the tool
``ase/ga/tools/get_all_candidates`` in the
same folder as the GA. This will create a trajectory file
``all_candidates.traj`` which includes all locally relaxed candidates
the GA has tried. This script can be run at the same time as the main
script is running. This is possible because the ase.db database
is being updated as the GA progresses.

Running the GA in Parallel
==========================
 
One of the great advantages of a GA is that many structures can be
relaxed in parallel. This GA implementation includes two classes which
facilitates running the GA in parallel. One class can be used for
running several single threaded optimizations simultaneously on the
same compute node, and the other class integrates the GA into the PBS
queuing system used at many high performance computer clusters.


Relaxations in Parallel on the Same Computer
--------------------------------------------

In order to relax several structures simultaneously on the same
computer a seperate script relaxing one structure needs to be
created. Continuing the example from above we therefore create a
script taking as input the filename of the structure to relax and
which as output saves a trajectory file with the locally optimized
structure. It is important that the relaxed structure is named as in
this script, since the parallel integration assumes this file naming
scheme. For the example described above this script could look like

.. literalinclude:: ga_basic_calc.py

The main script needs to initialize the parallel controller and then
the script needs to be changed the two places where structures are
relaxed. The changed main script now looks like

.. literalinclude:: ga_basic_parallel_main.py

Notice how the main script is not cluttered by the local optimization
logic and is therefore now also easier to read. ``n_simul`` controls
the number of simultaneous relaxations, and can of course also be set
to 1 effectively giving the same result as in the non parallel
situation.

The ``relax`` method on the ``ParallelLocalRun`` class only returns
control to the main script when there is an execution thread
available. In the above example the relax method immediately returns
control to the main script the first 4 times it is called, but the
fifth time control is first returned when one of the first four
relaxations have been completed.

Running the GA together with a queing system
============================================

The GA has been implemented with first principles structure
optimization in mind. When using for instance DFT calculations for the
local relaxations relaxing one structure can take many hours. For this
reason the GA has been made so that it can work together with queing
systems where each candidate is relaxed in a separate job. With this
in mind the main script of the GA can thus also be considered a
controller script which every time it is invoked gathers the current
population, checks with a queing system for the number of jobs
submitted, and submits new jobs. For a typical application the main
script can thus be invoked by a crontab once every hour.

To run the GA together with a queing system the user needs to specify
a function which takes as input a job name and the path to the
trajectory file that needs to be submitted (the ``jtg`` function in
the sample script below). From this the function generates a PBS job
file which is submitted to the queing system. The calculator script
specified in the jobfile needs to obey the same naming scheme as the
sample calculator script in the previous section. The sample
relaxation script given in the previous can be used as starting point
for a relaxation script.

Handling of the parallel logic is in this case in the main script. The
parameter n_simul given to the ``PBSQueueRun`` object determines how
many relaxations should be in the queuing system simultaneously. The
main script now looks the following:

.. literalinclude:: ga_basic_pbs_main.py
