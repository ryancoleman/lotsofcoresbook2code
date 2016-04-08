.. _dcdft_exercise:

==========================================================================
DeltaCodesDFT - Comparing Solid State DFT Codes, Basis Sets and Potentials
==========================================================================

The webpage https://molmod.ugent.be/deltacodesdft provides a method
for measuring the precision of a given calculation method
against a chosen reference method (computational
or experimental) for parameters of the equation of state
(see :ref:`aluminium_exercise`) of elementary solids.

When performing any benchmark calculations, especially involving
a large number of systems, it is important to be aware of the fact
that we, humans tend to do mistakes.

Therefore the motto of this exercise is taken from Karl Popper's
"All life is problem solving":
`the novelty in the scientific approach is that we actively seek to eliminate our attempted solutions <http://books.google.dk/books?id=W0jP04qn0uoC&pg=PA9&lpg=PA9&dq=%22The+novelty+in+the+scientific+approach+is+that+we+actively+seek+to+eliminate+our+attempted+solutions%22&source=bl&ots=pvgW-0uZUp&sig=Rr6pyMIoFa7Fq_RxiHkcdvjAgyo&hl=en&sa=X&ei=_kLSU_noJKWGywOLmYG4Aw&ved=0CCAQ6AEwAA#v=onepage&q=%22The%20novelty%20in%20the%20scientific%20approach%20is%20that%20we%20actively%20seek%20to%20eliminate%20our%20attempted%20solutions%22&f=false>`_.

In this exercise, in addition to the traditional, error prone
method of writing output files generated using separate scripts
on disk, we write the results into a database.
All calculations are performed using one script, and therefore
not only sharing the results with other researches is easy
(by granting the access to the database) but also the precise
method of performing the calculations should be shared (by presenting
the script). Please consult the introduction to the
:mod:`ase.db` module for details.

You can find out more about "reproducible science" with ASE in
the following talks:
`Emacs + org-mode + python in reproducible research
<http://www.youtube.com/watch?v=1-dUkyn_fZA>`_ or
`How Python & the iPython notebook can revamp quantum chemical reseach
<http://www.youtube.com/watch?v=WKoImDmYFQE>`_.

We will compare PBE numbers from GPAW with http://www.wien2k.at/ for K, Ca
and Ti. We use default PAW-datasets, a plane-wave cutoff of 340 eV, a
**k**-point density of 3.5 Å and a Fermi-Dirac distribution with a 0.1 eV
temperature. Copy this :download:`dcdft_gpaw.py` to a place in your file area:

.. literalinclude:: dcdft_gpaw.py

.. highlight:: bash

Read the script and try to understand it. Run the
script by typing::

  $ python dcdft_gpaw.py

It should take about 15 minutes to run the script.
Note that you can start several instances of the script
simultaneously in order to speed things up.

The script will generate ``.txt`` files and an SQLite3 database file.
Watch the progess as the calculations run::

    $ ase-db dcdft.db -c +x,time

Examine the equation of state (see :ref:`aluminium_exercise`)
using :command:`ase-gui`::

  $ ase-gui dcdft.db@name=Ca

.. note::

    The PBE reference values from https://molmod.ugent.be/deltacodesdft are:
     
    =======  ==================  =========
    element  `V` [Å\ `^3`/atom]  `B` [GPa]
    =======  ==================  =========
    K        73.68               3.6
    Ca       42.20               17.1
    Ti       17.39               112.2
    =======  ==================  =========

Extract the results from the database in order to calculate
the parameters of the equation of state:

.. literalinclude:: dcdft.db_raw.txt

and use the script available from https://molmod.ugent.be/deltacodesdft
to calculate the Delta factors.

* How well do the obtained values agree with the references?
  Do you think they can be further improved
  (hint: check out https://wiki.fysik.dtu.dk/gpaw/setups/dcdft.html)?
  Do you agree with Karl Popper?
