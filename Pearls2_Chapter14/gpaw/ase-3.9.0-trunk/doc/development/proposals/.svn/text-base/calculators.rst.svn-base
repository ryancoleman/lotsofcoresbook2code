.. _aep1:

=============================
Calculator interface proposal
=============================

All ASE calculators should behave similarly if there is no good reason
for them not to.  This should make it simpler for both users and developers.

This proposal tries to define how a good ASE calculator should behave.
The goal is to have ASE calculators:

* that share more code
* are more uniform to use
* are better tested
* are portable

Setting some standards is a good thing, but we should also be careful
not to set too strict rules that could limit each calculator to the
lowest common denominator.


Behavior
========

When a calculator calculates the energy, forces, stress tensor, total
magnetic moment, atomic magnetic moments or dipole moment, it should
store a copy of the system (atomic numbers, atomic positions, unit
cell and boundary conditions).  When asked again, it should return the
value already calculated if the system hasn't been changed.

If calculational parameters such as plane wave cutoff or XC-functional
has been changed the calculator should throw away old calculated
values.


Standards parameters
====================

The standard keywords that all calculators must use (if they make
sense) are: ``xc``, ``kpts``, ``smearing``, ``charge`` and ``nbands``.
Each calculator will have its own default values for these parameters
--- see recommendations below.  In addition, calculators will
typically have many other parameters.  The units are eV and Å.

Initial magnetic moments are taken from the :class:`~ase.atoms.Atoms`
object.

:xc:

  It is recommended that ``'LDA'`` and ``'PBE'`` are valid options.

:kpts:

  * ``(1,1,1)``: Gamma-point
  
  * ``(n1,n2,n3)``: Monkhorst-Pack grid
  
  * ``(n1,n2,n3,'gamma')``: Shifted Monkhorst-Pack grid that includes `\Gamma`
  
  * ``[(k11,k12,k13),(k21,k22,k23),...]``: Explicit list in units of the
    reciprocal lattice vectors
  
  * ``kpts=3.5``: `\vec k`-point density as in 3.5 `\vec k`-points per
    Å\ `^{-1}`

:smearing:

  The smearing parameter must be given as a tuple:

  * ``('Fermi-Dirac', width)``
  * ``('Gaussian', width)``
  * ``('Methfessel-Paxton', width, n)``, where `n` is the order (`n=0`
    is the same as ``'Gaussian'``)

  Lower-case strings are also allowed.  The ``width`` parameter used
  for the chosen smearing method is in eV units.

:charge:

  Charge of the system in units of `|e|` (``charge=1`` means one
  electron has been removed).


:nbands:

  Each band can be occupied by two electrons.

  
ABC calculator example
======================

The constructor will look like this::

  ABC(restart=None, ignore_bad_restart=False, label=None,
      atoms=None, **kwargs)

A calculator should be able to prefix all output files with a given
label or run the calculation in a directory with a specified name.
This is handled by the ``label`` argument.  There are three
possibilities:

* Name of a file containing all results of a calculation (possibly
  containing a directory).

* A prefix used for several files containing results.  The label may
  have both a directory part and a prefix part like ``'LDA/mol1'``.

* Name of a directory containing result files with fixed names.

Each calculator can decide what the default value is: ``None`` for no
output, ``'-'`` for standard output or something else.

If the ``restart`` argument is given, atomic configuration, input
parameters and results will be read from a previous calculation from
the file(s) pointed to by the ``restart`` argument.  It is an error if
those files don't exist and are corrupted.  This error can be ignored
bu using ``ignore_bad_restart=True``.

The ``atoms`` argument is discussed below.  All additional parameters
are given as keyword arguments.

Example:  Do a calculation with ABC calculator and write results to
:file:`si.abc`:

>>> atoms = ...
>>> atoms.calc = ABC(label='si.abc', xc='LDA', kpts=3.0)
>>> atoms.get_potential_energy()
-1.2

Read atoms with ABC calculator attaced from a previous calculation:

>>> atoms = ABC.read_atoms('si.abc')
>>> atoms.calc
<ABC-calculator>
>>> atoms.get_potential_energy()
-1.2

The ``ABC.read_atoms('si.abc')`` statement is equivalent to::

  ABC(restart='si.abc', label='si.abc').get_atoms()

If we do:

>>> atoms = ABC.read_atoms('si.abc')
>>> atoms.rattle()            # change positions and/or
>>> atoms.calc.set(xc='PBE')  # change a calculator-parameter
>>> atoms.get_potential_energy()
-0.7

then the :file:`si.abc` will be overwritten or maybe appended to.

An alternative way to connect atoms and calculator:

>>> atoms = ...
>>> calc = ABC(restart='si.abc', label='si.abc', atoms=atoms)
>>> atoms.get_potential_energy()
-0.7

This will automatically attach the calculator to the atoms and the
atoms will be updated form the file.  If you add
``ignore_bad_restart=True``, you will be able to use the same
script to do the initial calculation where :file:`si.abc` does not
exist and following calculations where atoms may have been moved
arround by an optimization algorithm.

The command used to start the ABC code can be given in an environment
variable called :envvar:`ASE_ABC_COMMAND` or as a ``command``
keyword.  The command can look like this::

  mpiexec abc PREFIX.input > PREFIX.output

or like this::

  ~/bin/start_abc.py PREFIX

The ``PREFIX`` strings will be substituted by the ``label`` keyword.


Implementation
==============

* Portability (Linux/Windows): ``os.system('Linux commands')`` not allowed.

* Common base class for all calculators: ``Calculator``.  Takes care
  of restart from file logic, handles setting of parameters and checks
  for state changes.

* A ``FileIOCalculator`` for the case where we need to:

  * write input file(s)
  * run Fortran/C/C++ code
  * read output file(s)

* Helper function to deal with ``kpts`` keyword.
