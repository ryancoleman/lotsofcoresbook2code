.. _manual:

======
Manual
======

GPAW calculations are controlled through scripts written in the
programming language Python_.  GPAW relies on the :ase:`Atomic
Simulation Environment <>` (ASE), which is a Python package that helps
us describe our atoms.  The ASE package also
handles molecular dynamics, analysis, visualization, geometry
optimization and more.  If you don't know anything about ASE, then it
might be a good idea to familiarize yourself with it before continuing
(at least read the :ase:`ASE introduction <intro.html>`).

Below, there will be Python code examples starting with ``>>>`` (and
``...`` for continuation lines).  It is a good idea to start the
Python interpreter and try some of the examples below.

.. _Python: http://www.python.org

The units used by the GPAW calculator correspond to the :ase:`ASE
conventions <ase/units.html>`, most importantly electron volts and
angstroms.

-----------------------
Doing a PAW calculation
-----------------------

To do a PAW calculation with the GPAW code, you need an ASE
:class:`~ase.atoms.Atoms` object and a :class:`~gpaw.aseinterface.GPAW`
calculator::

   _____________          ____________
  |             |        |            |
  | Atoms       |------->| GPAW       |
  |             |        |            |
  |_____________|        |____________|
       atoms                  calc

In Python code, it looks like this:

.. literalinclude:: h2.py

If the above code was executed, a calculation for a single `\rm{H}_2`
molecule would be started.  The calculation would be done using a
supercell of size :math:`6.0 \times 6.0 \times 6.0` Å with cluster
boundary conditions.  The parameters for the PAW calculation are:

* 2 electronic bands.
* Local density approximation (LDA)\ [#LDA]_ for the
  exchange-correlation functional.
* Spin-paired calculation.
* :math:`32 \times 32 \times 32` grid points.

The values of these parameters can be found in the text output file:
`h2.txt <../h2.txt>`_.

The calculator will try to make sensible choices for all parameters
that the user does not specify.  Specifying parameters can be done
like this:

>>> calc = GPAW(nbands=1,
...             xc='PBE',
...             gpts=(24, 24, 24))

Here, we want to use one electronic band, the Perdew, Burke, Ernzerhof
(PBE)\ [#PBE]_ exchange-correlation functional and 24 grid points in
each direction.


----------
Parameters
----------

The complete list of all possible parameters and their defaults is
shown below. A detailed description of the individual parameters is
given in the following sections.

=================  =========  ===================  ============================
keyword            type       default value        description
=================  =========  ===================  ============================
``mode``           ``str``    ``'fd'``             :ref:`manual_mode`
``nbands``         ``int``                         :ref:`manual_nbands`
``xc``             ``str``    ``'LDA'``            :ref:`manual_xc`
``kpts``           *seq*      `\Gamma`-point       :ref:`manual_kpts`
``spinpol``        ``bool``                        :ref:`manual_spinpol`
``gpts``           *seq*                           :ref:`manual_gpts`
``h``              ``float``  ``0.2``              :ref:`manual_h`
``symmetry``       ``dict``   ``{}``               :ref:`manual_symmetry`
``usesymm``        ``bool``   ``True``             :ref:`manual_symmetry`
``random``         ``bool``   ``False``            Use random numbers for
                                                   :ref:`manual_random`
``occupations``    occ. obj.                       :ref:`manual_occ`
``lmax``           ``int``    ``2``                Maximum angular momentum
                                                   for expansion of
                                                   :ref:`manual_lmax`
``charge``         ``float``  ``0``                Total :ref:`manual_charge`
                                                   of the system
``convergence``    ``dict``                        :ref:`manual_convergence`
``maxiter``        ``int``    ``333``              :ref:`manual_maxiter`
``txt``            ``str``,   ``'-'``              :ref:`manual_txt`
                   None, or   (``sys.stdout``)
                   file obj.
``parallel``       ``dict``                        :ref:`manual_parallel`
``stencils``       tuple      ``(3, 3)``           Number of neighbors for
                                                   :ref:`manual_stencils`
``mixer``          Object                          Pulay :ref:`manual_mixer`
                                                   scheme
``fixdensity``     ``bool``   ``False``            Use :ref:`manual_fixdensity`
``setups``         ``str``    ``'paw'``            :ref:`manual_setups`
                   or
                   ``dict``
``basis``          ``str``    ``{}``               Specification of
                   or                              :ref:`manual_basis`
                   ``dict``
``eigensolver``    ``str``    ``'rmm-diis'``       :ref:`manual_eigensolver`
``hund``           ``bool``   ``False``            :ref:`Use Hund's rule
                                                   <manual_hund>`
``external``       Object                          XXX Missing doc
``verbose``        ``int``    ``0``                :ref:`manual_verbose`
``poissonsolver``  Object                          Specification of
                                                   :ref:`Poisson solver 
                                                   <manual_poissonsolver>`
                                                   or :ref:`dipole correction
                                                   <manual_dipole_correction>`
``communicator``   Object                          :ref:`manual_communicator`
``idiotproof``     ``bool``   ``True``             Set to ``False`` to ignore 
                                                   setup fingerprint mismatch
                                                   (allows restart when the
                                                   original setup files are
                                                   not available)
=================  =========  ===================  ============================

*seq*: A sequence of three ``int``'s.


.. note:: 
   
   Parameters can be changed after the calculator has been constructed
   by using the :meth:`~gpaw.aseinterface.GPAW.set` method:

   >>> calc.set(txt='H2.txt', charge=1)

   This would send all output to a file named :file:`'H2.txt'`, and the
   calculation will be done with one electron removed.


Deprecated keywords (in favour of the ``parallel`` keyword) include:

=================  =========  ===================  ============================
keyword            type       default value        description
=================  =========  ===================  ============================
``parsize``        *seq*                           Parallel
                                                   :ref:`manual_parsize_domain`
``parsize_bands``  ``int``    ``1``                :ref:`manual_parsize_bands`
=================  =========  ===================  ============================


.. _manual_mode:

Finite-difference, plane-wave or LCAO mode
------------------------------------------

Finite-difference:
    The default mode (``mode='fd'``) is Finite Differece. This means that
    the wave functions will be expanded on a real space grid.
    
LCAO:
    Expand the wave functions in a basis-set constructed
    from atomic-like orbitals, in short LCAO (linear combination of atomic
    orbitals).  This is done by setting ``mode='lcao'``.

    See also the page on :ref:`lcao`.

Plane-waves:
    Expand the wave-functions in plane-waves.  Use ``mode='pw'`` if you want
    to use the default plane-wave cutoff of `E_{\text{cut}}=340` eV.  The
    plane-waves will be those with `|\bG+\bk|^2/2<E_{\text{cut}}`.  You
    can set another cutoff like this::
        
        from gpaw import GPAW, PW
        calc = GPAW(mode=PW(200))


Comparing FD, LCAO and PW modes
```````````````````````````````
    
Memory consumption:
    With LCAO, you have fewer degrees of freedom so memory usage is low.
    PW mode uses more memory and FD a lot more.
    
Speed:
    For small systems with many **k**-points, PW mode beats everything else.
    For larger systems LCAO will be most efficient.  Whereas PW beats FD for
    smallish systems, the opposite is true for very large systems where FD
    will parallelize better.
    
Absolute convergence:
    With LCAO, it can be hard to reach the complete basis set limit and get
    absolute convergence of energies, whereas with FD and PW mode it is
    quite easy to do by decreasing the grid spacing or increasing the
    plane-wave cutoff energy, respectively.
    
Eggbox errors:
    With LCAO and FD mode you will get a small eggbox error: you get a
    small periodic energy variation as you translate atoms and the period
    of the variation will be equal to the grid-spacing used.  GPAW's PW
    implementation doesn't have this problem.

Features:
    FD mode is the oldest and has most features.  Only PW mode can be used
    for calculating the stress-tensor and for response function calculations.
    
    
.. _manual_nbands:

Number of electronic bands
--------------------------

This parameter determines how many bands are included in the calculation for
each spin. For example, for spin-unpolarized system with 10 valence electrons
``nbands=5`` would include all the occupied states. In 10 valence electron
spin-polarized system with magnetic moment of 2 a minimum of ``nbands=6`` is
needed (6 occupied bands for spin-up, 4 occupied bands and 2 empty bands for
spin down).

The default number of electronic bands (``nbands``) is equal to the
number of atomic orbitals present in the atomic setups.  For systems
with the occupied states well separated from the unoccupied states,
one could use just the number of bands needed to hold the occupied
states.  For metals, more bands are needed.  Sometimes, adding more
unoccupied bands will improve convergence.

.. tip::
   
    ``nbands=0`` will give zero empty bands, and ``nbands=-n`` will
    give ``n`` empty bands.
   
.. tip::
    
    ``nbands='n%'`` will give ``n/100`` times the number of occupied bands.


.. _manual_xc:

Exchange-Correlation functional
-------------------------------

Some of the most commonly used exchange-correlation functionals
are listed below.

============  =========================== ===========================  ==========
``xc``        full libxc_ keyword         description                  reference 
============  =========================== ===========================  ==========
``'LDA'``     ``'LDA_X+LDA_C_PW'``        Local density approximation  [#LDA]_
``'PBE'``     ``'GGA_X_PBE+GGA_C_PBE'``   Perdew, Burke, Ernzerhof     [#PBE]_
``'revPBE'``  ``'GGA_X_PBE_R+GGA_C_PBE'`` revised PBE                  [#revPBE]_
``'RPBE'``    ``'GGA_X_RPBE+GGA_C_PBE'``  revised revPBE               [#RPBE]_
``'PBE0'``    ``'HYB_GGA_XC_PBEH'``       Known as PBE0                [#PBE0]_
``'B3LYP'``   ``'HYB_GGA_XC_B3LYP'``      B3LYP (as in Gaussian Inc.)  [#B3LYP]_
============  =========================== ===========================  ==========

``'LDA'`` is the default value.  The next three ones are of
generalized gradient approximation (GGA) type, and the last two are
`hybrid functionals <http://en.wikipedia.org/wiki/Hybrid_functional>`_.

For the list of all functionals available in GPAW see :ref:`overview_xc`.

GPAW uses the functionals from libxc_ by default.
Keywords are based on the names in the libxc :file:`'xc_funcs.h'` header file (the leading ``'XC_'`` should be removed from those names).
Valid keywords are strings or combinations of exchange and correlation string
joined by **+** (plus).
For example, "the" (most common) LDA approximation in chemistry
corresponds to ``'LDA_X+LDA_C_VWN'``.

Hybrid functionals (the feature is described at :ref:`exx`)
require the setups containing exx information to be generated.
Check available setups for the presence of exx information, for example::

     [~]$ zcat $GPAW_SETUP_PATH/O.PBE.gz | grep "<exact_exchange_"

and generate setups with missing exx information::

     [~]$ gpaw-setup --exact-exchange -f PBE H C

Currently all the hybrid functionals use the PBE setup as a *base* setup.

For more information about ``gpaw-setup`` see :ref:`generation_of_setups`.

Set the location of setups as decribed on :ref:`installationguide_setup_files`.

The details of the implementation of the exchange-correlation
are described on the :ref:`xc_functionals` page.

.. _libxc: http://www.tddft.org/programs/octopus/wiki/index.php/Libxc


.. _manual_kpts:

Brillouin-zone sampling
-----------------------

The default sampling of the Brillouin-zone is with only the
`\Gamma`-point.  This allows us to choose the wave functions to be
real.  Monkhorst-Pack sampling can be used if required: ``kpts=(N1,
N2, N3)``, where ``N1``, ``N2`` and ``N3`` are positive integers.
This will sample the Brillouin-zone with a regular grid of ``N1``
`\times` ``N2`` `\times` ``N3`` **k**-points.  See the
:func:`ase.dft.kpoints.monkhorst_pack` function for more details.

For more flexibility, you can use this syntax::
    
    kpts={'size': (4, 4, 4)}  # 4x4x4 Monkhorst-pack
    kpts={'size': (4, 4, 4), 'gamma': True}  # shifted 4x4x4 Monkhorst-pack

You can also specify the **k**-point density in units of points per
Å\ `^{-1}`::
    
    kpts={'density': 2.5}  # Monkhorst-Pack with a density of 2.5 points/Ang^-1
    kpts={'density': 2.5, 'even': True}  # round off to neares even number
    kpts={'density': 2.5, 'gamma': True}  # include gamma-point
    
The **k**-point density is calculated as:
    
.. math:: N \frac{a}{2\pi},

where `N` is then number of **k**-points and `a` is the length of the
unit-cell along the direction of the corresponding reciprocal lattice vector.
    
An arbitrary set of **k**-points can be specified, by giving a
sequence of k-point coordinates like this::

    kpts=[(0, 0, -0.25), (0, 0, 0), (0, 0, 0.25), (0, 0, 0.5)]

The **k**-point coordinates are given in scaled coordinates, relative
to the basis vectors of the reciprocal unit cell.

The above four **k**-points are equivalent to
``kpts={'size': (1, 1, 4), 'gamma': True}`` and to this:

>>> from ase.dft.kpoints import monkhorst_pack
>>> kpts = monkhorst_pack((1, 1, 4))
>>> kpts
array([[ 0.   ,  0.   , -0.375],
       [ 0.   ,  0.   , -0.125],
       [ 0.   ,  0.   ,  0.125],
       [ 0.   ,  0.   ,  0.375]])
>>> kpts+=(0,0,0.125)
>>> kpts
array([[ 0.  ,  0.  , -0.25],
       [ 0.  ,  0.  ,  0.  ],
       [ 0.  ,  0.  ,  0.25],
       [ 0.  ,  0.  ,  0.5 ]])


.. _manual_spinpol:

Spinpolarized calculation
-------------------------

If any of the atoms have magnetic moments, then the calculation will
be spin-polarized - otherwise, a spin-paired calculation is carried
out.  This behavior can be overruled with the ``spinpol`` keyword
(``spinpol=True``).


.. _manual_gpts:

Number of grid points
---------------------

The number of grid points to use for the grid representation of the
wave functions determines the quality of the calculation.  More
gridpoints (smaller grid spacing, *h*), gives better convergence of
the total energy.  For most elements, *h* should be 0.2 Å for
reasonable convergence of total energies.  If a ``n1`` `\times` ``n2``
`\times` ``n3`` grid is desired, use ``gpts=(n1, n2, n3)``, where
``n1``, ``n2`` and ``n3`` are positive ``int``'s all divisible by four.
Alternatively, one can use something like ``h=0.25``, and the program
will try to choose a number of grid points that gives approximately
a grid-point density of `1/h^3`.  For more details, see :ref:`grids`.

If you are more used to think in terms of plane waves; a conversion
formula between plane wave energy cutoffs and realspace grid spacings
have been provided by Briggs *et. al* PRB **54**, 14362 (1996).  The
conversion can be done like this::

  >>> from gpaw.utilities.tools import cutoff2gridspacing, gridspacing2cutoff
  >>> from ase.units import Rydberg
  >>> h = cutoff2gridspacing(50 * Rydberg)


.. _manual_h:

Grid spacing
------------

The parameter ``h`` specifies the grid spacing in Å that has to be
used for the realspace representation of the smooth wave
functions. Note, that this grid spacing in most cases is approximate
as it has to fit to the unit cell (see :ref:`manual_gpts` above).

In case you want to specify ``h`` exactly you have to choose the unit
cell accordingly. This can be achieved by::

  from gpaw.cluster import *

  d = 0.74
  a = 6.0
  atoms = Cluster('H2', positions=[(0, 0, 0), (0, 0, d)])
  # set the amount of vacuum at least to 4 Å 
  # and ensure a grid spacing of h=0.2
  atoms.minimal_box(4., h=.2)


.. _manual_symmetry:

Use of symmetry
---------------

The default behavior is to use all point-group symmetries and time-reversal
symmetry to reduce the **k**-points to only those in the irreducible part of
the Brillouin-zone.  Moving the atoms so that a symmetry is broken will
cause an error.  This can be avoided by using::
    
    symmetry={'point_group': False}

This will reduce the number of applied symmetries to just the time-reversal
symmetry (implying that the Hamiltonian is invariant under **k** -> -**k**).
For some purposes you might want to have no symmetry reduction of the
**k**-points at all (debugging, transport calculations, band-structure
calculations, ...). This can be achieved by specifying::

    symmetry={'point_group': False, 'time_reversal': False}
    
or simply ``symmetry='off'`` which is a short-hand notation for the same
thing.

For full control, here are all the available keys of the ``symmetry``
dictionary:

=================  ========  ===============================
key                default   description
=================  ========  ===============================
``point_group``    ``True``  Use point-group symmetries
``time_reversal``  ``True``  Use time-reversal symmetry
``symmorphic``     ``True``  Use only symmorphic symmetries
``tolerance``      ``1e-7``  Relative tolerance
=================  ========  ===============================

.. note::
    
    If you are using version 0.10 or earlier, you can use
    ``usesymm=False`` to turn off all point-group symmetries and
    ``usesymm=None`` to turn off also time-reversal symmetry.
    

.. _manual_random:

Wave function initialization
----------------------------

By default, a linear combination of atomic orbitals is used as initial
guess for the wave functions. If the user wants to calculate more bands
than there are precalculated atomic orbitals, random numbers will be
used for the remaining bands.


.. _manual_occ:

Occupation numbers
------------------

The smearing of the occupation numbers is controled like this::

  from gpaw import GPAW, FermiDirac
  calc = GPAW(..., occupations=FermiDirac(width), ...)

The distribution looks like this (width = `k_B T`):

.. math::  f(E) = \frac{1}{1 + \exp[E / (k_B T)]}

For calculations with periodic boundary conditions, the default value
is 0.1 eV and the total energies are extrapolated to *T* = 0 Kelvin.
For a molecule (no periodic boundaries) the default value is ``width=0``,
which gives integer occupation numbers.

For a spin-polarized calculation, one can fix the magnetic moment at
the initial value using ``FermiDirac(width, fixmagmom=True)``.


.. _manual_lmax:

Compensation charges
--------------------

The compensation charges are expanded with correct multipoles up to
and including `\ell=\ell_{max}`.  Default value: ``lmax=2``.


.. _manual_charge:

Charge
------

The default is charge neutral.  The systems total charge may be set in
units of the negative electron charge (i.e. ``charge=-1`` means one
electron more than the neutral).


.. _manual_convergence:

Accuracy of the self-consistency cycle
--------------------------------------

The ``convergence`` keyword is used to set the convergence criteria.
The default value is this Python dictionary::

  {'energy': 0.0005,  # eV / electron
   'density': 1.0e-4,
   'eigenstates': 4.0e-8,  # eV^2 / electron
   'bands': 'occupied',
   'forces': None} # eV / Ang Max 

In words:

* The energy change (last 3 iterations) should be less than 0.5 meV
  per valence electron.

* The change in density (integrated absolute value of density change) 
  should be less than 0.001 electrons per valence electron.

* The integrated value of the square of the residuals of the Kohn-Sham
  equations should be less than :math:`4.0 \times 10^{-8}
  \mathrm{eV}^2` per valence electron (FD mode only).

* The maximum change in the magnitude of the vector representing the
  difference in forces for each atom.  Setting this to None disables 
  this functionality, saving computational time and memory usage.

The individual criteria can be changed by giving only the specific
entry of dictionary e.g. ``convergence={'energy': 0.0001}`` would set
the convergence criteria of energy to 0.1 meV while other criteria
remain in their default values.

As the total energy and charge density depend only on the occupied
states, unoccupied states do not contribute to the convergence
criteria.  However, with the ``bands`` set to ``'all'``, it is
possible to force convergence also for the unoccupied states.  One can
also use ``{'bands': 200}`` to converge the lowest 200 bands. One can
also write ``{'bands': -10}`` to converge all bands except the last
10. It is often hard to converge the last few bands in a calculation.


.. _manual_maxiter:

Maximum number of SCF-iterations
--------------------------------

The calculation will stop with an error if convergence is not reached
in ``maxiter`` self-consistent iterations (defaults to 120).


.. _manual_txt:

Where to send text output
-------------------------

The ``txt`` keyword defaults to the string ``'-'``, which means
standard output.  One can also give a ``file`` object (anything with a
``write`` method will do).  If a string (different from ``'-'``) is
passed to the ``txt`` keyword, a file with that name will be opened
and used for output.  Use ``txt=None`` to disable all text output.


.. _manual_mixer:

Density mixing
--------------

Three parameters determine how GPAW does Pulay mixing of the
densities:

* ``beta``: linear mixing coefficient
* ``nmaxold``: number of old densities to mix
* ``weight``: when measuring the change from input to output density,
  long wavelength changes are weighted ``weight`` times higher than
  short wavelength changes

For small molecules, the best choice is to use
``mixer=Mixer(beta=0.25, nmaxold=3, weight=1.0)``, which is what GPAW
will choose if the system has zero-boundary conditions.

If your system is a big molecule or a cluster, it is an advantage to
use something like ``mixer=Mixer(beta=0.1, nmaxold=5, weight=50.0)``,
which is also what GPAW will choose if the system has periodic
boundary conditions in one or more directions.

In spin-polarized calculations using Fermi-distribution
occupations one has to use :class:`~gpaw.mixer.MixerSum` instead of
:class:`~gpaw.mixer.Mixer`.

See also the documentation on :ref:`density mixing <densitymix>`.

.. _manual_fixdensity:

Fixed density
-------------

When calculating band structures or when adding unoccupied states to
calculation (and wanting to converge them) it is often useful to use existing
density without updating it. By using ``fixdensity=True`` the initial density 
(e.g. one read from .gpw/.hdf5 or existing from previous calculation) is used
throughout the SCF-cycle (so called Harris calculation).



.. _manual_setups:

Type of setup to use
--------------------

The ``setups`` keyword is used to specify the name(s) of the setup files
used in the calulation.

For a given element ``E``, setup name ``NAME``, and xc-functional
'XC', GPAW looks for the file :file:`E.NAME.XC` or :file:`E.NAME.XC.gz`
(in that order) in the setup locations
(see :ref:`installationguide_setup_files`).
Unless ``NAME='paw'``, in which case it will simply look for
:file:`E.XC` (or :file:`E.XC.gz`).
The ``setups`` keyword can be either a single string, or a dictionary.

If specified as a string, the given name is used for all atoms.  If
specified as a dictionary, each keys can be either a chemical symbol
or an atom number. The values state the individual setup names.

The special key ``None`` can be used to specify the default setup
name. Thus ``setups={None: 'paw'}`` is equivalent to ``setups='paw'``
which is the GPAW default.

As an example, the latest PAW setup of Na includes also the 6 semicore p states
in the valence, in order to use non-default setup with only the 1 s electron in valence (:file:`Na.1.XC.gz`) one can specify ``setups={'Na': '1'}``

There exist three special names, that if used, does not specify a file name:

* ``'ae'`` is used for specifying all-electron mode for an
  atom. I.e. no PAW or pseudo potential is used.
* ``'hgh'`` is used to specify a norm-conserving Hartwigsen-Goedecker-Hutter 
  pseudopotential (no file necessary).  Some elements have better 
  semicore pseudopotentials.  To use those, specify ``'hgh.sc'`` 
  for the elements or atoms in question.
* ``'ghost'`` is used to indicated a *ghost* atom in LCAO mode, 
  see :ref:`ghost-atoms`. 

If a dictionary contains both chemical element specifications *and*
atomic number specifications, the latter is dominant.

An example::

  setups={None: 'soft', 'Li': 'hard', 5: 'ghost', 'H': 'ae'}
  
Indicates that the files named 'hard' should be used for lithium
atoms, an all-electron potential is used for hydrogen atoms, atom
number 5 is a ghost atom (even if it is a Li or H atom), and for all
other atoms the files named 'soft' will be used.


.. _manual_basis:

Atomic basis set
----------------

The ``basis`` keyword can be used to specify the basis set which
should be used in LCAO mode, which also affects the LCAO
initialization in FD mode.

In FD mode, the initial guess for the density / wave functions is
determined by solving the Kohn-Sham equations in the LCAO basis.

The ``basis`` keyword can be either a string or a dictionary.  If
``basis`` is a string, GPAW will look for a file named
:file:`{symbol}.{basis}.basis` in
the setup locations
(see :ref:`installationguide_setup_files`), where
:file:`{symbol}` is taken as the chemical symbol from the ``Atoms``
object.
The first found file is used.
If ``basis`` is a dictionary, the basis set can be specified
differently for each atomic species by using the atomic symbol as
a key, or for individual atoms by using an ``int`` as a key.  In the
latter case the integer corresponds to the index of that atom in the
``Atoms`` object.  As an example, ``basis={'H': 'sz', 'C': 'dz', 7:
'dzp'}`` will use the ``sz`` basis for hydrogen atoms, the ``dz``
basis for carbon, and the ``dzp`` for whichever atom is number 7 in
the ``Atoms`` object.

.. note::

    If you want to use only the ``sz`` basis functinons from a ``dzp``
    basis set, then you can use this syntax: ``basis='sz(dzp)'``.
    This will read the basis functions for, say hydrogen, from the
    ``H.dzp.basis`` file.

The value ``None`` (default) implies that the pseudo partial waves
from the setup are used as a basis. This basis is always available;
choosing anything else requires the existence of the corresponding
basis set file in setup locations
(see :ref:`installationguide_setup_files`).

For details on the LCAO mode and generation of basis set files; see
the :ref:`LCAO <lcao>` documentation.


.. _manual_eigensolver:

Eigensolver
-----------

The default solver for iterative diagonalization of the Kohn-Sham
Hamiltonian is RMM-DIIS (Residual minimization method - direct
inversion in iterative subspace) which seems to perform well in most
cases. However, some times more efficient/stable convergence can be
obtained with a different eigensolver. Especially, when calculating
many unoccupied states RMM-DIIS might not be optimal. Further
available options in FD mode are conjugate gradient method
(``eigensolver='cg'``) and a simple Davidson method
(``eigensolver='dav'``). From the alternatives, conjugate gradient
seems to perform better in general.

More control can be obtained by using directly the eigensolver objects::

  from gpaw.eigensolvers import CG
  calc = GPAW(eigensolver=CG(niter=5, rtol=0.20))

Here, ``niter`` specifies the maximum number of conjugate gradient iterations
for each band (within a single SCF step), and if the relative change 
in residual is less than ``rtol``, the iteration for the band is not continued.

LCAO mode has its own eigensolver, which directly diagonalizes the
Hamiltonian matrix instead of using an iterative method.


.. _manual_poissonsolver:

Poisson solver
--------------

The *poissonsolver* keyword is used to specify a Poisson solver class
or enable dipole correction.

The default Poisson solver uses a multigrid Jacobian method.  Use this
keyword to specify a different method.  This example corresponds to
the default Poisson solver::

  from gpaw import GPAW, PoissonSolver
  calc = GPAW(poissonsolver=PoissonSolver(nn=3, relax='J', eps=2e-10))

The first argument is the stencil, see :ref:`manual_stencils`.  Second
argument is the method, either ``'J'`` (Jacobian) or ``'GS'``
(Gauss-Seidel).  The Gauss-Seidel method requires half as many
iterations to solve the Poisson equation, but involves more
communication.  The Gauss-Seidel implementation also depends slightly
on the domain decomposition used.

The last argument, ``eps``, is the convergence criterion.

.. note:: 

  The Poisson solver is rarely a performance bottleneck, but it can
  sometimes perform poorly depending on the grid layout.  This is mostly 
  important in LCAO calculations, but can be good to know in general.
  See the LCAO notes on
  :ref:`Poisson performance <poisson_performance>`.

.. _manual_dipole_correction:

The *poissonsolver* keyword can also be used to specify that a dipole
correction should be applied along a given axis.  The system should be
non-periodic in that direction but periodic in the two other
directions.

::
  
  from gpaw import GPAW
  from gpaw.poisson import PoissonSolver
  from gpaw.dipole_correction import DipoleCorrection
  
  poissonsolver = PoissonSolver()
  correction = DipoleCorrection(poissonsolver, 2) # 2 == z-axis
  calc = GPAW(poissonsolver=correction)

Without dipole correction, the potential will approach 0 at all
non-periodic boundaries.  With dipole correction, there will be a
potential difference across the system depending on the size of the
dipole moment.


.. _manual_stencils:

Finite-difference stencils
--------------------------

GPAW uses finite-difference stencils for the Laplacian in the
Kohn-Sham and Poisson equations.  You can set the range of the stencil
(number of neighbor grid points) used for the Poisson equation like
this::

    from gpaw import GPAW, PoissonSolver
    calc = GPAW(poissonsolver=PoissonSolver(nn=n))

This will give an accuracy of `O(h^{2n})`, where ``n`` must be between
1 and 6.  The default value for version 0.6 is ``n='M'`` which is a
special Mehrstellen stencil - this will be changed to ``n=3`` in
version 0.7.

With the ``stencils=(a, b)`` keyword, you can set the accuracy of the
stencil used for the Kohn-Sham equation to `O(h^{2a})`.  The ``b``
parameter (between 1 and 4) controls the accuracy of the
interpolation of the density from the coarse grid to the fine grid.
Default values are ``stencils=(3, 3)``.


.. _manual_hund:

Using Hund's rule for guessing initial magnetic moments
-------------------------------------------------------

The ``hund`` keyword can be used for single atoms only. If set to
``True``, the calculation will become spinpolarized, and the initial
ocupations, and magnetic moment of the atom will be set to the value
required by Hund's rule.  You may further wish to specify that the
total magnetic moment be fixed, by passing e.g.
``occupations=FermiDirac(0.0, fixmagmom=True)``.
Any user specified magnetic moment is
ignored. Default is False.

.. _manual_verbose:

Output verbosity
----------------

By default, only a limited number of information is printed out for each SCF
step. It is possible to obtain more information (e.g. for investigating 
convergen problems in more detail) by ``verbose=1`` keyword.

.. _manual_communicator:

Communicator object
-------------------

By specifying a communicator object, it is possible to use only a subset of
processes for the calculator when calculating e.g. different atomic images 
in parallel. See :ref:`different_calculations_in parallel` for more details.

.. _manual_parallel_calculations:

---------------------
Parallel calculations
---------------------

Information about running parallel calculations can be found on the
:ref:`parallel_runs` page.


.. _zero_energy:

--------------
Total Energies
--------------

The GPAW code calculates energies relative to the energy of separated
reference atoms, where each atom is in a spin-paired, neutral, and
spherically symmetric state - the state that was used to generate the
setup.  For a calculation of a molecule, the energy will be minus the
atomization energy and for a solid, the resulting energy is minus the
cohesive energy.  So, if you ever get positive energies from your
calculations, your system is in an unstable state!

.. note::
   You don't get the true atomization/cohesive energy.  The true
   number is always lower, because most atoms have a spin-polarized
   and non-spherical symmetric ground state, with an energy that is
   lower than that of the spin-paired, and spherically symmetric
   reference atom.


------------------------
Restarting a calculation
------------------------

The state of a calculation can be saved to a file like this:

>>> calc.write('H2.gpw')

The file :file:`H2.gpw` is a binary file containing
wave functions, densities, positions and everything else (also the
parameters characterizing the PAW calculator used for the
calculation).

If you want to restart the `\rm{H}_2` calculation in another Python session
at a later time, this can be done as follows:

>>> from gpaw import *
>>> atoms, calc = restart('H2.gpw')
>>> print atoms.get_potential_energy()

Everything will be just as before we wrote the :file:`H2.gpw` file.
Often, one wants to restart the calculation with one or two parameters
changed slightly.  This is very simple to do.  Suppose you want to
change the number of grid points:

>>> atoms, calc = restart('H2.gpw', gpts=(20, 20, 20))
>>> print atoms.get_potential_energy()

.. tip::
   There is an alternative way to do this, that can be handy sometimes:

   >>> atoms, calc = restart('H2.gpw')
   >>> calc.set(gpts=(20, 20, 20))
   >>> print atoms.get_potential_energy()


More details can be found on the :ref:`restart_files` page.

---------------------------------------
Customizing behaviour through observers
---------------------------------------

An *observer* function can be *attached* to the calculator so that it
will be executed every *N* iterations during a calculation.  The below
example saves a differently named restart file every 5 iterations::

  calc = GPAW(...)

  occasionally = 5

  class OccasionalWriter:
      def __init__(self):
          self.iter = 0

      def write(self):
          calc.write('filename.%03d.gpw' % self.iter)
          self.iter += occasionally

  calc.attach(OccasionalWriter().write, occasionally)

See also :meth:`~gpaw.aseinterface.GPAW.attach`.


----------------------
Command line arguments
----------------------

The behaviour of GPAW can be controlled with some command line
arguments. The arguments for GPAW should be specified after the
python-script, i.e.::

    python script.py [options]

The possible command line arguments are:

===============================  =============================================
argument                         description
===============================  =============================================
``--trace``
``--debug``                      
                                 Run in debug-mode, e.g. check
                                 consistency of arrays passed to c-extensions
``--dry-run[=nprocs]``           
                                 Print out the computational
                                 parameters and estimate memory usage, 
                                 do not perform actual calculation. 
                                 If ``nprocs`` is specified, print also which 
                                 parallelization settings would be employed
                                 when run on ``nprocs`` processors.
``--memory-estimate-depth[=n]``  
                                 Print out an itemized memory estimate by
                                 stepping recursively through the object
                                 hierarchy of the calculator. If ``n`` is
                                 specified, print a summary for depths
                                 greater than the specified value.
                                 Default: ``n=2``
``--domain-decomposition=comp``  
                                 Specify the domain decomposition with
                                 ``comp`` as a positive integer or, for
                                 greater control, a tuple of three integers.
                                 Allowed values are equivalent to those of
                                 the ``domain`` argument in the
                                 :ref:`parallel <manual_parallel>` keyword,
                                 with tuples specified as ``nx,ny,nz``.
                                 See :ref:`manual_parsize_domain` for details.
``--state-parallelization=nbg``
                                 Specify the parallelization over Kohn-Sham
                                 orbitals with ``nbg`` as a positive integer.
                                 Allowed values are equivalent to those of
                                 the ``band`` argument in the
                                 :ref:`parallel <manual_parallel>` keyword.
                                 See :ref:`manual_parsize_bands` for details.
``--sl_...=m,n,mb``
                                 Specify ScaLAPACK / BLACS parameters for
                                 diagonalization (``--sl_default``),
                                 inverse Cholesky factorization
                                 (``--sl_inverse_cholesky``) and LCAO general
                                 diagonalization (``--sl_lcao``) respectively.
                                 Use ``--sl_default`` to specify all three
                                 of the above at once or just the default
                                 value for those not explicitly given.
                                 Allowed values are equivalent to those of
                                 the four ``sl_...`` arguments in the 
                                 :ref:`parallel <manual_parallel>` keyword.
                                 Requires GPAW to be built with ScaLAPACK.
``--gpaw a=1,b=2.3,...``         
                                 Extra parameters for development work:
                                 
                                 >>> from gpaw import extra_parameters
                                 >>> print extra_parameters
                                 {'a': 1, 'b': 2.3}
===============================  =============================================


.. [#LDA]    J. P. Perdew and Y. Wang,
             Accurate and simple analytic representation of the
             electron-gas correlation energy
             *Phys. Rev. B* **45**, 13244-13249 (1992)
.. [#PBE]    J. P. Perdew, K. Burke, and M. Ernzerhof,
             Generalized Gradient Approximation Made Simple,
             *Phys. Rev. Lett.* **77**, 3865 (1996)
.. [#revPBE] Y. Zhang and W. Yang,
             Comment on "Generalized Gradient Approximation Made Simple",
             *Phys. Rev. Lett.* **80**, 890 (1998)
.. [#RPBE]   B. Hammer, L. B. Hansen and J. K. Nørskov,
             Improved adsorption energetics within density-functional
             theory using revised Perdew-Burke-Ernzerhof functionals,
             *Phys. Rev. B* **59**, 7413 (1999)
.. [#PBE0]   C. Adamo and V. Barone,
             *J. Chem. Phys.* **110** 6158-6170 (1999)
             Toward reliable density functional methods without adjustable
             parameters: The PBE0 model
.. [#B3LYP]  P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and M.J. Frisch,
             *J. Phys. Chem.* **98** 11623-11627 (1994)
             Ab-Initio Calculation of Vibrational Absorption and Circular-Dichroism
             Spectra Using Density-Functional Force-Fields
