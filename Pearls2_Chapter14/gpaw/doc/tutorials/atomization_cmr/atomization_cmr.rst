.. _atomization_energy_cmr:

===================================================
Calculation of atomization energies with aid of CMR
===================================================

The goal of this tutorial is to illustrate the use of
Computational Materials Repository (CMR_) in calculation
of atomization energies. It is assumed you are familiar
with :ref:`atomization_energy`.

.. _CMR:
    https://wiki.fysik.dtu.dk/cmr/

The full script :svn:`gpaw/test/cmrtest/Li2_atomize.py`
is divided into several parts:

1. the preliminary part contains important settings relevant for CMR:

  - setting a **unique** project identifier:

    .. literalinclude:: ../../../gpaw/test/cmrtest/Li2_atomize.py
       :start-after: project
       :end-before: vacuum

    .. note::

       You will access your results by ``project_id``,
       so don't forget to set it!

  - defining the parameters of the calculation to be stored in CMR:

    .. literalinclude:: ../../../gpaw/test/cmrtest/Li2_atomize.py
       :start-after: template
       :end-before: }

2. ``calculate`` part, starting with::

     if calculate:

   setups the Li2 molecule and the Li atom, and performs ``LDA``
   (see :ref:`manual_xc`) total energy calculations.
   The calculations are saved to the ``Li2.gpw`` and ``Li.gpw``
   files (see :ref:`restart_files`)::

     calc.write(formula)

   Additionally, the required parameters of the calculations
   are saved into the corresponding cmr files::

     write(cmrfile, system, cmr_params=cmr_params)

3. ``recalculate`` part, starting with::

     if recalculate:

   performs calculation of the ``PBE`` total energy on ``LDA`` orbitals.
   The orbitals are read from the :file:`.gpw` files saved in the ``calculate`` step.
   The results are appended to the corresponding cmr files:

   .. literalinclude:: ../../../gpaw/test/cmrtest/Li2_atomize.py
      :start-after: exists
      :end-before: del

Please perform the steps up to (and including) the ``recalculate`` part.

.. note::

   In order to do so set the control variables accordingly::

     calculate = True
     recalculate = True
     analyse_from_dir = False # analyse local cmr files

     upload_to_db = False  # upload cmr files to the database
     analyse_from_db = False # analyse database

     create_group = False # group calculations belonging to a given reaction

     clean = False

You can now calculate atomization energy of the Li2 molecule by opening
the text output files and extracting the results.

In the next part the atomization energy is calculated based on
the cmr files saved in the current directory.

4. the ``analyse_from_dir`` part, starting with::

     if analyse_from_dir:

   performs the analysis.

   .. note::

      Set the control variables accordingly::

        calculate = False
        recalculate = False
        analyse_from_dir = True # analyse local cmr files

        upload_to_db = False  # upload cmr files to the database
        analyse_from_db = False # analyse database

        create_group = True # group calculations belonging to a given reaction

        clean = False

   In this part the contents of all cmr files in the current directory is read,
   and restricted to our ``project_id``:

   .. literalinclude:: ../../../gpaw/test/cmrtest/Li2_atomize.py
      :start-after: DirectoryReader
      :end-before: rank

   The ``LDA`` and ``PBE`` (on ``LDA`` orbitals) atomization energies are
   calculated with, respectively:

   .. literalinclude:: ../../../gpaw/test/cmrtest/Li2_atomize.py
      :start-after: (ea)
      :end-before: print

   and a group is created in order to connect the result of the calculation
   to the cmr files the calculation is based on.

   The result is::

     atomization energy [eV] LDA = 1.10
     atomization energy [eV] PBE = 0.99

   to be compared against :ref:`molecule_tests`.

5. all the results can be uploaded to a database.
   This is performed in the ``upload_to_db`` part, starting with::

     if upload_to_db:

   After waiting few minutes (time it takes to upload the results),
   one can calculate the atomization energies by querying the database directly.

   .. note::

      As an idication that the calculations are performed on the results
      from the database, please remove all the cmr files from the current directory!

6. analysis is performed by the ``analyse_from_db`` part.

