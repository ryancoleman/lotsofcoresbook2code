==============================
Initialization and I/O changes
==============================

This is a proposal for some changes that will solve various issues with 
the maintainability and stability of the I/O code amongst other things.

.. contents::

Rationale
=========

Presently the gpw I/O is handled centrally by the module 
gpaw/io/__init__.py.  If someone makes changes in setup.py or 
density.py, the I/O may break due to these "non-local correlations" (we 
in particular, being physicists, should know to appreciate locality), or 
it may misbehave in subtle ways for certain cases 
(TDDFT/LCAO/non-gamma-point/etc.).

Most of this trouble can be avoided entirely by requiring that objects 
should know how to read and write themselves.  Thus, responsibility for 
what to write (and how to read it back!) is delegated to various objects 
as per the usual 'object oriented' way.

A different but sort of related issue: The output.py module writes lots 
of things to the log file, and those things would be better off 'writing 
themselves'.  There are several bugs here waiting to be fixed: if the 
Poisson solver was created by hand with a non-standard stencil, then the 
wrong stencil is written.  Scalapack/BLACS information is sometimes 
wrong (depending on the way it was specified).  No information on the 
stridedness of the band descriptor is written (and thus parallelization 
info is incomplete).  There are probably other issues.


Object hierarki
===============

.. image:: ../bigpicture.png
   :target: ../../bigpicture.svg

So all the objects above may implement functions to read and write their 
own parameters.  They could also implement functions to read/write 
human-readable information to log files (which is highly desirable).

On a somewhat deeper level, we could formalize the tree hierarchy 
between the major GPAW objects and define a mandatory interface for all 
major objects to implement (read/write, memory estimate, 
initialize/set_positions/allocate -- these procesure *all* involve 
traversing the same tree hierarchy).  This might make it easier for new 
programmers to learn the structure of GPAW (a well-known problem).

Example of what an object could look like:

.. literalinclude:: density.py


The PAW calculator object
=========================

The following base class (rough sketch) should propably be moved to
ASE, so that all ASE-calculators can use it:

.. literalinclude:: ase.py

It should be possible to create the PAW calculator without knowing the
atoms and also from a restart file - and both ways must be cheap.

.. literalinclude:: paw.py

Open questions
--------------

The above pseudo code is not the final answer - it is an attempt to
make the further discussions more concrete.  There are several things
to think about:

* What should the ASECalculator look like?
* How much should/can ``__init__()`` for the different objects do?


Reading and writing
===================

We should name things like described here: http://www.etsf.eu/fileformats

(is there an ETSF XC library?)

Things we need to write:

* Version and units.

* Atoms: atomic numbers, positions, initial magnetic moments, tags,
  boundary conditions, unit cell, charges.

* Setups: lmax, fingerprints, setup types.

* Basis set.

* Hamiltonian: Poisson solver, XC functional, effective
  pseudopotential, non-local part of hamiltonian.

* Density: Charge, density convergence criterion, density error,
  atomic density matrices, interpolation order, pseudoelectron
  density, multipole moments of compensation charges.

* Mixer.

* Occupations: fixed magnetic moment flag, smearing type, width, Fermi
  level, occupation numbers.

* Symmmetry: Symmetry matrices, atom maps.

* Parameters that the result should not depend on: hund, random,
  maxiter, eigensolver?, parallel (domain, band, stridebands, scalapack).

* SCF: energy convergence criterion, energy error.

* Eigensolver: eigenstates convergence criterion, number of bands to
  converge, eigenstate error(s)?

* Calculated stuff: Ekin, Epot, Ebar, Eext, Exc, S, forces, potential
  energy, magnetic moments, dipole moment.

* Brillouin zone sampling: BZ, IBZ, weights, maps, symmetries.

* Projections.

* Pseudo wave functions.

* Eigenvalues.

What do we do with these: fixdensity, mode, Kohn Sham stencil, h, charge?

What do we need in order to better support:

* DSCF
* DFPT
* response functions
* GW
* TDDFT
* NEGF transport
* LCAO
* plane waves
* other things?

How should reading and writing work?

* Should it be like pickle where the class name is written (example:
  gpaw.mixer.MixerSum).  The reader will then create that object and
  read into it?

* What methods should a reader have?  ``get_object('name of
  object')``, ``get_array('name of array')``, ``get_atoms()``,
  ``get_grid_descriptor()``, ...

* We need backwards compatibility.

Also, there should be well defined interface so that it is easy to use
different backends (.gpw, hdf5, ...). I think that the current
io-interface ('__set_item__', 'dimension', 'add', 'fill', ...) is
quite well defined, and if new IO scheme requires
additions/modifications, they should be such that different backends
can easily support them.


Some thoughts about different backends:
--------------------------------------- 

If/when each object writes/reads itself, some sort of hierachical file
format would be convenient.  I am not that familiar with
tarfile-interface used for .gpw files, but I think that it can support
hierachical structure (folders and files). Also, HDF5 supports
hierachical structure ("hierachical data format"), basic structure of
HDF5 file is groups and datasets.  Other formats that one could think
of are MPI-IO and netcdf, but that they do not really support
hierarchical structure. Drawback of MPI-IO is also that the files are
not necessarily portable (although it should be possible to ensure
portability with the price of more expensive IO).


Here is a prototype implementation of a hierachical reader/writer
framework: :trac:`doc/devel/proposals/rw.py`.


Parallel IO
===========

For large calculations it will be more or less compulsory
to perform IO in parallel.  Even though a filesystem would not support
parallel IO (meaning that read/write are not faster than in serial
case), memory requirements can prohibit collecting the data into
single process. As an example, in large calculation with e.g. 200**3
grid, collecting density into single process requires 8 * 400**3 ~ 500
MB.

Some backends supporting parallel IO are MPI-IO, parallel-netcdf, and
HDF5, and there are existing python interfaces to MPI-IO (mpi4py) and
HDF5 (h5py and pytables). GPAW can already use h5py without parallel
capabilities. Enabling parallel IO with h5py is quite simple as it
requires adding only two simple functions to GPAW.

At some point, we should start using mpi4py with GPAW.

Backends
========

Tarfile
-------

Relatively simple, portable and no external dependencies, but:

* no parallel IO, single process has to collect data
* no direct access with external software (visualization etc.)

HDF5
----

Portable, can perform parallel IO and external software can access the
file directly, but:

* additional dependencies (at least HDF5 library, a python interface
  could in principle be included in GPAW)
* porting to more exotic architectures (Cray, Blue Gene) can be tricky?

Directory
---------

A bit like an extraxted tarfile.  Different cpu's could write
different states.  When the writing is done, one can tar the directory
to get a standard gpw file.  The tarfile format would have to be
modifyed so that one can read pseudo wave functions from several
files.
