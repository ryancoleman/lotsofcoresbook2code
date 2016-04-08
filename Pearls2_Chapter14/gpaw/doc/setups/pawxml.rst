.. _pawxml:

=========================================
XML specification for atomic PAW datasets
=========================================

------------
Introduction
------------

This page contains information about the PAW-XML data format for the
atomic datasets necessary for doing Projector Augmented-Wave
calculations \ [#Blo94]_.  We use the term *dataset* instead of
*pseudo potential* because the PAW method is not a pseudopotential method.

An example XML file for nitrogen PAW dataset using LDA can be seen
here: `N.LDA <../N.LDA>`_.

.. note::
   Hartree atomic units are used in the XML file (`\hbar = m = e = 1`).


-----------------------
What defines a dataset?
-----------------------

The following quantities defines a minimum PAW dataset (the notation
from Ref. [#Blo03]_ is used here):

============================  ======================================
Quantity                      Description
============================  ======================================
`Z`                           atomic number
`E_\text{XC}[n]`              exchange-correlation functional
`E^\text{kin}_c`              kinetic energy of the core electrons
`g_{\ell m}(\mathbf{r})`      shape function for compensation charge
`n_c(r)`                      all-electron core density
`\tilde{n}_c(r)`              pseudo electron core density
`\tilde{n}_v(r)`              pseudo electron valence density
`\bar{v}(r)`                  zero potential
`\phi_i(\mathbf{r})`          all-electron partial waves
`\tilde{\phi}_i(\mathbf{r})`  pseudo partial waves
`\tilde{p}_i(\mathbf{r})`     projector functions
`\Delta E^\text{kin}_{ij}`    kinetic energy differences
============================  ======================================

The following quantities can be optionally provided:

============================  ===============================================
Quantity                      Description
============================  ===============================================
`r_{PAW}`                     Radius of the PAW augmentation region (max. of matching radii)
`v_H[\tilde{n}_{Zc}](r)`      Kresse-Joubert local ionic pseudopotential
`\tilde{Q}_{ij}(\mathbf{r})`  State-dependent shape function for compensation charge
`\tau_c(r)`                   Core kinetic energy density
`\tilde{\tau}_c(r)`           Pseudo core kinetic energy density
`X^{\text{core-core}}`        Core-core contribution to exact exchange
`X_{ij}^{\text{core-val}}`    Core-valence exact-exchange correction matrix
============================  ===============================================


-----------------------------
Specification of the elements
-----------------------------

An element looks like this::

  <name> ... </name>

or for an empty element::

  <name/>

.. tip::
   An XML-tutorial can be found here_

   .. _here: http://www.w3schools.com/xml/default.asp


----------
The header
----------

The first two lines should look like this::

  <?xml version="1.0"?>
  <paw_dataset version="0.7">

The first line must be present in all XML files.  Everything else is put
inside an element with name ``paw_dataset``, and this element has an
attribute called ``version``.  We are currently at version 0.7.


---------
A comment
---------

It is recommended to put a comment giving the units and a link to this
web page::

  <!-- Nitrogen dataset for the Projector Augmented Wave method. -->
  <!-- Units: Hartree and Bohr radii.                            -->
  <!-- http://www.where.org/paw_dataset.html                     -->


--------------------
The ``atom`` element
--------------------

::

    <atom symbol="N" Z="7" core="2" valence="5"/>

The ``atom`` element has attributes ``symbol``, ``Z``, ``core`` and
``valence`` (chemical symbol, atomic number, number of core electrons and
number of valence electrons).


--------------------
Exchange-correlation
--------------------

The ``xc_functional`` element defines the exchange-correlation
functional used for generating the dataset. It has the two attributes ``type`` and ``name``.

The ``type`` attribute can be ``LDA``, ``GGA``, ``MGGA`` or ``HYB``.

The ``name`` attribute designates the exchange-correlation functional and
can be specified in the following ways:
 
- Taking the names from the LibXC_ library. The correlation and exchange names are stripped
  from their ``XC_`` part and combined with a ``+``-sign.  Here is an
  example for an LDA functional::
    
  <xc_functional type="LDA", name="LDA_X+LDA_C_PW"/>

  and this is what PBE will look like::

  <xc_functional type="GGA", name="GGA_X_PBE+GGA_C_PBE"/>

- Using one of the following pre-defined aliases:

  =========  ==========  ===============================  ===================================================================
  ``type``    ``name``   ``LibXC equivalent``             ``Reference``
  =========  ==========  ===============================  ===================================================================
  ``LDA``    ``PW``      ``LDA_X+LDA_C_PW``               ``LDA exchange; Perdew, Wang, PRB 45, 13244 (1992)``
  ``GGA``    ``PW91``    ``GGA_X_PW91+GGA_C_PW91``        ``Perdew et al PRB 46, 6671 (1992)``
  ``GGA``    ``PBE``     ``GGA_X_PBE+GGA_C_PBE``          ``Perdew, Burke, Ernzerhof, PRL 77, 3865 (1996)``
  ``GGA``    ``RPBE``    ``GGA_X_RPBE+GGA_C_PBE``         ``Hammer, Hansen, Nørskov, PRB 59, 7413 (1999)``
  ``GGA``    ``revPBE``  ``GGA_X_PBE_R+GGA_C_PBE``        ``Zhang, Yang, PRL 80, 890 (1998)``
  ``GGA``    ``PBEsol``  ``GGA_X_PBE_SOL+GGA_C_PBE_SOL``  ``Perdew et al, PRL 100, 136406 (2008)``
  ``GGA``    ``AM05``    ``GGA_X_AM05+GGA_C_AM05``        ``Armiento, Mattsson, PRB 72, 085108 (2005)``
  ``GGA``    ``BLYP``    ``GGA_X_B88+GGA_C_LYP``          ``Becke, PRA 38, 3098 (1988); Lee, Yang, Parr, PRB 37, 785 (1988)``
  =========  ==========  ===============================  ===================================================================

  Examples::
  
    <xc_functional type="LDA", name="PW"/>

  ::

  <xc_functional type="GGA", name="PBE"/>

.. _LibXC: http://www.tddft.org/programs/octopus/wiki/index.php/
           Libxc:manual#Available_functionals


---------
Generator
---------

::

  <generator type="scalar-relativistic" name="MyGenerator-2.0">
    Frozen core: [He]
  </generator>


This element contains *character data* describing in words how the
dataset was generated.  The ``type`` attribute must be one of:
``non-relativistic``, ``scalar-relativistic`` or ``relativistic``.


--------
Energies
--------

::

  <ae_energy kinetic="53.777460" xc="-6.127751"
             electrostatic="-101.690410" total="-54.040701"/>
  <core_energy kinetic="43.529213"/>

The kinetic energy of the core electrons,
`E^\text{kin}_c`, is used in the PAW method.  The other
energies are convenient to have for testing purposes and can also be
useful for checking the quality of the underlying atomic calculation.


--------------
Valence states
--------------

::

  <valence_states>
    <state n="2" l="0" f="2"  rc="1.10" e="-0.6766" id="N-2s"/>
    <state n="2" l="1" f="3"  rc="1.10" e="-0.2660" id="N-2p"/>
    <state       l="0"        rc="1.10" e=" 0.3234" id="N-s1"/>
    <state       l="1"        rc="1.10" e=" 0.7340" id="N-p1"/>
    <state       l="2"        rc="1.10" e=" 0.0000" id="N-d1"/>
  </valence_states>

The ``valence_states`` element contains several ``state`` elements, defined by a unique ``id``
as well as ``l`` and ``n`` quantum numbers. For each of them it is also required to provide
the energy ``e``, the occupation ``f``
and the matching radius of the partial waves ``rc``.

The number of ``state`` elements determines the size of the partial wave basis.
It is equal to the number of `radial functions`__
(radial parts of the `\phi_i`, `\tilde{\phi}_i` and `\tilde{p}_i`)
and is noted `n_{waves}` in the rest of this document.

For this dataset, the first two lines describe bound eigenstates with
occupation numbers and principal quantum numbers.  Notice, that the
three additional unbound states should have no ``f`` and ``n``
attributes.  In this way, we know that only the first two bound states
(with ``f`` and ``n`` attributes) should be used for constructing an
initial guess for the wave functions.

__ `Radial functions`_


------------
Radial grids
------------

There can be one or more definitions of radial grids.

Example::

  <radial_grid eq="r=d*i" d="0.1" istart="0" iend="9" id="g1">
    <values>
      0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
    </values>
    <derivatives>
      0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
    </derivatives>
  </radial_grid>
    
This defines one radial grid as :math:`r_i = di` where `i` runs from 0 to 9.  Inside the ``<radial_grid>`` element we have the
10 values of `r_i` followed by the 10 values of the derivatives
`dr_i/di`.

All functions (densities, potentials, ...) that use this grid are given as 10 numbers defining
the radial part of the function.  The radial part of the function must
be multiplied by a spherical harmonics:
`f_{\ell m}(\mathbf{r}) = f_\ell(r) Y_{\ell m}(\theta, \phi)`.

Each radial grid has a unique id::

  <radial_grid eq="r=d*i" d="0.01" istart="0" iend="99" id="lin">
  <radial_grid eq="r=a*exp(d*i)" a="1.056e-4" d="0.05" istart="0" iend="249" id="log">

and each numerical function must refer to one of these ids::

  <function grid="lin">
    ... ... ...
  </function>

In this example, the ``function`` element should contain 100 numbers
(`i = 0, ..., 99`).  Each number must be separated by a ``<newline>``
character or by one or more ``<tab>``'s or ``<space>``'s (no commas).
For numbers with scientific notation, use this format: ``1.23456e-5``
or ``1.23456E-5`` and not ``1.23456D-5``.

A program can read the values for `r_i` and `dr_i/di` from the file or
evaluate them from the ``eq`` and associated parameter attributes.
There are currently six types of radial grids:

=====================  ========================
``eq``                 parameters
=====================  ========================
``r=d*i``              ``d``
``r=a*exp(d*i)``       ``a`` and ``d``
``r=a*(exp(d*i)-1)``   ``a`` and ``d``
``r=a*i/(1-b*i)``      ``a`` and ``b``
``r=a*i/(n-i)``        ``a`` and ``n``
``r=(i/n+a)^5/a-a^4``  ``a`` and ``n``
=====================  ========================

The ``istart`` and ``iend`` attributes indicating the range of `i`
should always be present.


Although it is possible to define as radial grids as desired, it is recommended
to minimize the number of grids in the dataset.


------------------------------------------
Shape function for the compensation charge
------------------------------------------

The general formulation of the compensation charge uses an expansion over the partial
waves *ij* and the spherical harmonics:

.. math::

  \sum_{\ell m} C_{\ell m \ell_i m_i \ell_j m_j} \hat{Q}^{\ell}_{i j}(r) Y_{\ell m}(\theta, \phi),


where :math:`C_{\ell m \ell_i m_i \ell_j m_j}` is a *Gaunt coefficient*.

The standard expression \ [#Blo94]_ for the *shape function* :math:`\hat{Q}^{\ell}_{i j}(\mathbf{r})`
is a product of the multipole moment :math:`Q^{\ell}_{i j}` and a shape function :math:`g_\ell(r)`:

.. math::

  \hat{Q}^{\ell}_{i j}(r) = Q^{\ell}_{i j} g_\ell(r),

Several formulations [#Hol01]_ [#Blo94]_ define
`g_\ell(r) \propto r^\ell k(r)`, where `k(r)` is an `\ell`-independent
shape function:

==========  ===================  =========================================
``type``    parameters           `k(r)`
==========  ===================  =========================================
``gauss``   ``rc``               `\exp(-(r/r_c)^2)`
``sinc``    ``rc``               `[\sin(\pi r/r_c)/(\pi r/r_c)]^2`
``exp``     ``rc`` and ``lamb``  `\exp(-(r/r_c)^\lambda)`
==========  ===================  =========================================

Example::
    
    <shape_function type="gauss" rc="3.478505426185e-01">

Another formulation [#Kre99]_ defines directly `g_\ell(r)`:

==========  ==========  ===============================================
``type``    parameters  `g_\ell(r)`
==========  ==========  ===============================================
``bessel``  ``rc``      `\sum_{i=1}^2 \alpha_i^\ell j_\ell(q_i^\ell r)`
==========  ==========  ===============================================

For ``bessel`` the four parameters (`\alpha_1^\ell`, `q_1^\ell`,
`\alpha_2^\ell` and `q_2^\ell`) must be determined from ``rc`` for each
value of `\ell` as described in [#Kre99]_.

Example::
    
    <shape_function type="bessel" rc="3.478505426185e-01">
 

There is also a more general formulation where :math:`\hat{Q}^{\ell}_{i j}(r)` is given in
a numerical form. Several *shape functions* can be set (with the ``<shape_function>`` tag),
depending on `\ell` and/or combinations of partial waves (specified using the optional
``state1`` and ``state2`` attributes).
See for instance section II.C of [#Laa93]_.

Example 1, defining numerically :math:`g_\ell(r)`
in :math:`\hat{Q}^{\ell}_{i j}(r)=Q^{\ell}_{i j} g_\ell(r)`::
    
    <shape_function type="numeric" l=0 grid="g1">
        ... ... ...
    </shape_function>


Example 2, defining directly :math:`\hat{Q}^{\ell}_{i j}(r)`
for states *i=* ``N-2s`` and *j=* ``N-2p``, and *l=0*::
    
    <shape_function type="numeric" l=0 state1="N-2s" state2="N-2p" grid="g1">
        ... ... ...
    </shape_function>


----------------
Radial functions
----------------

Continuing, we have now reached the *all-electron* (resp. *pseudo core*,
*pseudo valence*) density::

  <ae_core_density grid="g1">
     6.801207147443e+02 6.801207147443e+02 6.665042896724e+02
     ... ...
  </ae_core_density>
  <pseudo_core_density rc="1.1" grid="g1">
     ...
  </pseudo_core_density>
  <pseudo_valence_density rc="1.1" grid="g1">
     ...
  </pseudo_valence_density>

The numbers inside the ``ae_core_density`` (resp. ``pseudo_core_density``, ``pseudo_valence_density``)
element defines the radial part of `n_c(\mathbf{r})` (resp. `\tilde{n}_c(\mathbf{r})`,
`\tilde{n}_v(\mathbf{r})`).
The radial part must be multiplied by `Y_{00} = (4\pi)^{-1/2}` to get the full density.
(`Y_{00}n_c(\mathbf{r})` should integrate to the number of core electrons).
The *pseudo core density* and the *pseudo valence* density are defined similarly and also
have a ``rc`` attribute specifying the matching radius.
 

The ``ae_partial_wave``, ``pseudo_partial_wave`` and
``projector_function`` elements contain the radial parts of the
`\phi_i(\mathbf{r})`, `\tilde{\phi}_i(\mathbf{r})` and
`\tilde{p}_i(\mathbf{r})` functions for the ``state``\ s listed in
the ``valence_states`` element above (five states in the nitrogen
example).  All functions must have an attribute ``state="..."``
referring to one of the states listed in the ``valence_states``
element::

  <ae_partial_wave state="N-2s" grid="g1">
    -8.178800366898029e+00 -8.178246914143839e+00 -8.177654917302689e+00
    ... ...
  </ae_partial_wave>
  <pseudo_partial_wave state="N-2s" grid="g1">
    ...
  </pseudo_partial_wave>
  <projector_function state="N-2s" grid="g1">
    ...
  </projector_function>
  <ae_partial_wave state="N-2p" grid="g1">
    ...
  </ae_partial_wave>
  ...
  ...

Remember that the radial part of these functions must be multiplied by a spherical harmonics:
`\phi_i(\mathbf{r}) = \phi_i(r) Y_{\ell_i m_i}(\theta, \phi)`.


--------------------------
Zero potential
--------------------------

The zero potential, `\bar{v}` (see section VI.D of [#Blo94]_) is defined similarly to the
densities; the radial part must be multiplied by `Y_{00} = (4\pi)^{-1/2}` to get the full
potential. The ``zero_potential`` element has a ``rc`` attribute specifying the cut-off
radius of `\bar{v}(\mathbf{r})`::
 
  <zero_potential rc="1.1" grid="g1">
     ...
  </zero_potential>


------------------------------
The Kresse-Joubert formulation
------------------------------

The Kresse-Joubert formulation of the PAW method\ [#Kre99]_ is very
similar to the original formulation of Blöchl\ [#Blo94]_.
However, the Kresse-Joubert formulation does not use `\bar{v}`
directly, but indirectly through the local ionic pseudopotential,
`v_H[\tilde{n}_{Zc}]`.  Therefore, the following
transformation is necessary:

.. math::

  v_H[\tilde{n}_{Zc}] = v_H[\tilde{n}_c +
  (N_c - Z - \tilde{N}_c) g_{00} Y_{00}] + \bar{v} +
  v_{xc}[\tilde{n}_v + \tilde{n}_c] -
  v_{xc}[\tilde{n}_v + \tilde{n}_c +
         (N_v - \tilde{N}_v - \tilde{N}_c) g_{00} Y_{00}]

where `N_c` is the number of core electrons, `N_v` is the number of
valence electrons, `\tilde{N}_c` is the number of electrons contained
in the pseudo core density and `\tilde{N}_v` is the number of
electrons contained in the pseudo valence density.
The Hartree potential from the density `n` is defined as:

.. math::

   v_H[n](r_1) = 4\pi \int_0^\infty r_2^2 dr_2 \frac{n(r_2)}{r_>},

where `r_>` is the larger of `r_1` and `r_2`.

.. note::
   In the Kresse-Joubert formulation, the symbol `\tilde{n}` is used
   for what we here call `\tilde{n}_v` and in the Blöchl formulation,
   we have `\tilde{n} = \tilde{n}_c + \tilde{n}_v`.

It is also possible to add an element
``kresse_joubert_local_ionic_pseudopotential`` that contains the
`v_H[\tilde{n}_{Zc}](r)` function directly, so that no conversion is
necessary::

  <kresse_joubert_local_ionic_pseudopotential rc="1.3" grid="log">
     ...
  </kresse_joubert_local_ionic_pseudopotential>

The ``kresse_joubert_local_ionic_pseudopotential`` element has a ``rc`` attribute
specifying the matching radius. This matching radius corresponds to the maximum
of all the matching radii used in the formalism.


--------------------------
Kinetic energy differences
--------------------------

::

    <kinetic_energy_differences>
       1.744042161013e+00 0.000000000000e+00 2.730637956456e+00
       ...
    <kinetic_energy_differences>

This element contains the symmetric `\Delta E^\text{kin}_{ij}` matrix:

.. math::

  \Delta E^\text{kin}_{ij} = \langle \phi_i | \hat{T} | \phi_j \rangle
  - \langle \tilde{\phi}_i | \hat{T} | \tilde{\phi}_j \rangle

| where `\hat{T}` is the kinetic energy operator used by the generator.
| With `n_{waves}` valence states (see `n_{waves}` `definition`__), we have a `n_{waves} \times n_{waves}` matrix listed as `n_{waves}^2` numbers.

__ `Valence states`_


--------
Meta-GGA
--------

Datasets for use with MGGA functionals must also include information on the
*core kinetic energy density* and *pseudo core kinetic energy density* ;
the latters are defined with these two elements::
    
    <ae_core_kinetic_energy_density grid="g1">
      ... ... ...
    </ae_core_kinetic_energy_density>
    <pseudo_core_kinetic_energy_density rc="1.1" grid="g1">
      ... ... ...
    </pseudo_core_kinetic_energy_density>

These densities are defined similarly to the core and valence densities (see above).
The ``pseudo_core_kinetic_energy_density`` element has a ``rc`` attribute specifying its
matching radius.


------------------------
Exact exchange integrals
------------------------

The core-core contribution to the exact exchange energy
`X^{\text{core-core}}` and the symmetric core-valence
PAW-correction matrix `X_{ij}^{\text{core-valence}}` are given as:

.. math::
    
    X^{\text{core-core}} = -\frac{1}{4}\sum_{cc'} \iint d\mathbf{r} d\mathbf{r}'
    \frac{\phi_c(\mathbf{r})\phi_{c'}(\mathbf{r}) \phi_c(\mathbf{r}')\phi_{c'}(\mathbf{r}')}
    {|\mathbf{r}-\mathbf{r}'|}

.. math::

    X_{ij}^{\text{core-valence}} = -\frac{1}{2}\sum_c \iint d\mathbf{r} d\mathbf{r}'
    \frac{\phi_i(\mathbf{r})\phi_c(\mathbf{r}) \phi_j(\mathbf{r}')\phi_c(\mathbf{r}')}
    {|\mathbf{r}-\mathbf{r}'|}

The `X_{ij}^{\text{core-valence}}` coefficients depend only on pairs of the radial
basis functions `\phi_i(r)` and can be evaluated by summing over radial
integrals times **3-j** symbols according to:

.. math::

    X_{ij}^{\text{core-valence}} =
    -\delta_{l_i l_j} \delta_{m_i m_j} \sum_{c L} \frac{N_c}{2}
    {\begin{pmatrix}l_c & L & l_i \\ 0 & 0 & 0\end{pmatrix}}^2
    \int r^2 dr \int {r'}^2 d{r'}
    \frac{r^{L}_{<}}{r^{L+1}_{>}}
    \phi_i(r) \phi_c(r) \phi_j(r') \phi_c(r')

| where
| `N_{c}` is the number of core electrons corresponding to `l_{c}` (`N_c=2l_c+1`),
| `r_>` (resp. `r_<`) is the larger (resp. smaller) of `r` and `r'`.


`X^{\text{core-core}}` can be specified in the ``core`` attribute of the
``<exact_exchange>`` element.

 
With `n_{waves}` valence states (see `n_{waves}` `definition`__),
`X_{ij}^{\text{core-valence}}` is a `n_{waves} \times n_{waves}` matrix.
It can be specified as `n_{waves}^2` numbers inside the ``<exact_exchange>`` element::
    
    <exact_exchange core="...">
      ... ... ...
    </exact_exchange>

__ `Valence states`_


-----------------
Optional elements
-----------------

::

   <paw_radius rc="2.3456781234">

Although not necessary, it may be helpful to provide the following item(s) in the dataset:

 - Radius of the PAW augmentation region ``paw_radius``
   
   This radius defines the region (around the atom) outside which all pseudo quantities
   are equal to the all-electron ones.
   It is equal to the maximum of all the cut-off and matching radii.
   Note that -- for better lisibility -- the ``paw_radius`` element should be
   provided in the header of the file.


------------------
End of the dataset
------------------

::

  </paw_dataset>


-----------------------
How to use the datasets
-----------------------

Most likely, the radial functions will be needed on some other type of
radial grid than the one used in the dataset.  The idea is that one
should read in the radial functions and then transform them to the
radial grids used by the specific implementation.  After the
transformation, some sort of normalization may be necessary.


-----------------------------
Plotting the radial functions
-----------------------------

The first 10-20 lines of the XML-datasets, should be pretty much human
readable, and should give an overview of what kind of dataset it is and
how it was generated.  The remaining part of the file contains
numerical data for all the radial functions.  To get an overview of
these functions, you can extract that data with the
:trac:`~doc/setups/pawxml.py` program and then pass it on to your
favorite plotting tool.

.. note::
   The ``pawxml.py`` program is very primitive and is only included in
   order to demonstrates how to parse XML using SAX
   from a Python program.  Parsing XML from Fortran or C code with
   SAX should be similar.

Usage:

It works like this::

  $ pawxml.py [options] dataset[.gz]

Options:

==================================  =======================================
``--version``                       Show program's version number and exit.
``-h, --help``                      Show this help message and exit.
``-x <name>, --extract=<name>``     Function to extract.
``-s<channel>, --state=<channel>``  Select valence state.
``-l, --list``                      List valence states
==================================  =======================================

Examples::

  [~]$ pawxml.py -x pseudo_core_density N.LDA | xmgrace -
  [~]$ pawxml.py -x ae_partial_wave -s N2p N.LDA > N.ae.2p
  [~]$ pawxml.py -x pseudo_partial_wave -s N2p N.LDA > N.ps.2p
  [~]$ xmgrace N.??.2p


----------
References
----------

.. [#Blo94]  P. E. Blöchl,
             Projector augmented-wave method,
             *Phys. Rev. B* **50**, 17953-19979 (1994)
.. [#Kre99]  G. Kresse and D. Joubert,
             Form ultrasoft pseudopotentials to the projector
             augmented-wave method,
             *Phys. Rev. B* **59**, 1758-1775 (1999)
.. [#Hol01]  N. A. W. Holzwarth, A. R. Tackett, and G. E. Matthews,
             A Projector Augmented Wave (PAW) code for electronics
             structure calculations: Part I *atompaw* for generating
             atom-centered functions,
             *Computer Physics Communications* **135**, 329-347 (2001)
.. [#Blo03]  P. E. Blöchl, C. J. Forst and J. Schimpl,
             Projector augmented wave method: Ab initio molecular
             dynamics with full wave functions,
             *Bulletin of Materials Science* **26**, 33-41 (2003)
.. [#Laa93]  K. Laasonen, A. Pasquarello, R. Car, C. Lee and D. Vanderbilt,
             Car-Parrinello molecular dynamics with Vanderbilt
             ultrasoft pseudopotentials,
             *Phys. Rev. B* **47**, 10142-10153 (1993)
