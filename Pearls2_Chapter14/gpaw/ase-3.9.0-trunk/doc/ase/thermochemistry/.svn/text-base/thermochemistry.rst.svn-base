.. module:: ase.thermochemistry
   :synopsis: Thermochemistry module

===============
Thermochemistry
===============

ASE contains a :mod:`ase.thermochemistry` module that lets the user derive
commonly desired thermodynamic quantities of molecules and crystalline solids
from ASE output and some user-specified parameters. Three cases are currently
handled by this module: the ideal-gas limit (in which translational and
rotational degrees of freedom are taken into account), the harmonic limit
(generally used for adsorbates, in which all degrees of freedom are treated
harmonically), and a crystalline solid model (in which a lattice of N atoms
is treated as a system of 3N independent harmonic oscillators). The first two
cases rely on good vibrational energies being fed to the calculators, which
can be calculated with the :mod:`ase.vibrations` module. Likewise, the
crystalline solid model depends on an accurate phonon density of states; this
is readily calculated using the :mod:`ase.phonons` module.


Ideal-gas limit
===============

The thermodynamic quantities of ideal gases are calculated by assuming that
all spatial degrees of freedom are independent and separable into
translational, rotational, and vibrational degrees of freedom. The
:class:`~ase.thermochemistry.IdealGasThermo` class supports calculation of
enthalpy (:math:`H`), entropy (:math:`S`), and Gibbs free energy (:math:`G`),
and has the interface listed below.

.. autoclass:: IdealGasThermo
   :members:

       
Example
-------

The :class:`IdealGasThermo` class would generally be called after an energy
optimization and a vibrational analysis. The user needs to supply certain
parameters if the entropy or free energy are desired, such as the geometry
and symmetry number. An example on the nitrogen molecule is:

.. literalinclude:: nitrogen.py


This will give the thermodynamic summary output:

.. literalinclude:: nitrogen.txt


Harmonic limit
==============

In the harmonic limit, all degrees of freedom are treated harmonically. The
:class:`HarmonicThermo` class supports the calculation of internal energy,
entropy, and Gibbs free energy. This class uses all of the energies given to
it in the vib_energies list; this is a list as can be generated with the
.get_energies() method of :class:`ase.vibrations.Vibrations`, but the user
should take care that all of these energies are real (non-imaginary). The
class :class:`HarmonicThermo` has the interface described below.

.. autoclass:: HarmonicThermo
   :members:


Crystals
========

In this model a crystalline solid is treated as a periodic system of
independent harmonic oscillators. The :class:`CrystalThermo` class supports
the calculation of internal energy (:math:`U`), entropy (:math:`S`) and
Helmholtz free energy (:math:`F`), and has the interface listed below.

.. autoclass:: CrystalThermo
   :members:

       
Example
-------

The :class:`CrystalThermo` class will generally be called after an energy
optimization and a phonon vibrational analysis of the crystal. An example for
bulk gold is:

.. literalinclude:: gold.py

This will give the thermodynamic summary output:

.. literalinclude:: gold.txt


Background
==========

**Ideal gas.** The conversion of electronic structure calculations to
thermodynamic properties in the ideal-gas limit is well documented; see, for
example, Chapter 10 of Cramer, 2004. The key equations used in the
:class:`IdealGasThermo` class are summarized here.

   C.J. Cramer. *Essentials of Computational Chemistry*, Second Edition.
   Wiley, 2004.

The ideal-gas enthalpy is calculated from extrapolation of the energy at 0 K
to the relevant temperature (for an ideal gas, the enthalpy is not a function
of pressure):

.. math ::
   H(T) = E_\text{elec} + E_\text{ZPE} + \int_0^\text{T} C_P \, \text{d}T

where the first two terms are the electronic energy and the zero-point energy, and the integral is over the constant-pressure heat capacity. The heat capacity is separable into translational, rotational, vibrational, and electronic parts (plus a term of :math:`k_\text{B}` to switch from constant-volume to constant-pressure):

.. math ::
   C_P = k_\text{B} + C_{V\text{,trans}} + C_{V\text{,rot}} + C_{V\text{,vib}} + C_{V\text{,elec}}

The translational heat capacity is 3/2 :math:`k_\text{B}` for a 3-dimensional
gas. The rotational heat capacity is 0 for a monatomic species,
:math:`k_\text{B}` for a linear molecule, and 3/2 :math:`k_\text{B}` for a
nonlinear molecule. In this module, the electronic component of the heat
capacity is assumed to be 0. The vibrational heat capacity contains
:math:`3N-6` degrees of freedom for nonlinear molecules and :math:`3N-5`
degrees of freedom for linear molecules (where :math:`N` is the number of
atoms). The integrated form of the vibrational heat capacity is:

.. math ::
   \int_0^T C_{V,\text{vib}} \text{d}T = \sum_i^\text{vib DOF}
   \frac{\epsilon_i}{e^{\epsilon_i / k_\text{B} T} - 1 }

where :math:`\epsilon_i` are the energies associated with the vibrational
frequencies, :math:`\epsilon_i = h \omega_i`.

The ideal gas entropy can be calculated as a function of temperature and
pressure as:

.. math ::
   S(T,P) &= S(T,P^\circ) - k_\text{B} \ln \frac{P}{P^\circ} \\
          &= S_\text{trans} + S_\text{rot} + S_\text{elec} + S_\text{vib} - k_\text{B} \ln \frac{P}{P^\circ}

where the translational, rotational, electronic, and vibrational components
are calculated as below. (Note that the translational component also includes
components from the Stirling approximation, and that the vibrational degrees
of freedom are enumerated the same as in the above.)

.. math ::
   S_\text{trans} = k_\text{B} \left\{ \ln \left[ \left(
   \frac{2 \pi M k_\text{B} T}{h^2} \right)^{3/2}
   \frac{k_\text{B} T}{P^\circ} \right] + \frac{5}{2} \right\}

.. math ::
   S_\text{rot} = \left\{  \begin{array}{ll}
   0 & \text{, if monatomic} \\
   k_\text{B} \left[ \ln \left( \frac{8\pi^2 I k_\text{B}T}{\sigma h^2}\right) + 1 \right] & \text{, if linear} \\
   k_\text{B} \left\{ \ln \left[ \frac{\sqrt{\pi I_\text{A} I_\text{B} I_\text{C}}}{\sigma} \left(\frac{8\pi^2 k_\text{B} T}{h^2}\right)^{3/2}\right] + \frac{3}{2} \right\} & \text{, if nonlinear} \\
   \end{array}
   \right.
   
.. math ::
   S_\text{vib} = k_\text{B} \sum_i^\text{vib DOF}
   \left[ \frac{\epsilon_i}{k_\text{B}T\left(e^{\epsilon_i/k_\text{B}T}-1\right)} - \ln \left( 1 - e^{-\epsilon_i/k_\text{B}T} \right)\right]

.. math ::
   S_\text{elec} = k_\text{B} \ln \left[
   2 \times \left(\text{spin multiplicity}\right) + 1\right]

:math:`I_\text{A}` through :math:`I_\text{C}` are the three principle moments
of inertia for a non-linear molecule. :math:`I` is the degenerate moment of
inertia for a linear molecule. :math:`\sigma` is the symmetry number of the
molecule.

The ideal-gas Gibbs free energy is then just calculated from the combination
of the enthalpy and entropy:

.. math ::
   G(T,P) = H(T) - T\, S(T,P)

**Harmonic limit.** The conversion of electronic structure calculation
information into thermodynamic properties is less established for adsorbates.
However, the simplest approach often taken is to treat all :math:`3N` degrees
of freedom of the adsorbate harmonically since the adsorbate often has no
real translational or rotational degrees of freedom. This is the approach
implemented in the :class:`HarmonicThermo` class. Thus,
the internal energy and entropy of the adsorbate are calculated as

.. math ::
   U(T) = E_\text{elec} + E_\text{ZPE} + \sum_i^\text{harm DOF} \frac{\epsilon_i}{e^{\epsilon_i / k_\text{B} T} - 1 }

.. math ::
   S = k_\text{B} \sum_i^\text{harm DOF}
   \left[ \frac{\epsilon_i}{k_\text{B}T\left(e^{\epsilon_i/k_\text{B}T}-1\right)} - \ln \left( 1 - e^{-\epsilon_i/k_\text{B}T} \right)\right]

and the Gibbs free energy is calculated as

.. math ::
   G(T) = U(T) - T\, S(T)

In this case, the number of harmonic energies (:math:`\epsilon_i`) used in
the summation is generally :math:`3N`, where :math:`N` is the number of atoms
in the adsorbate.

**Crystalline solid**

The derivation of the partition function for a crystalline solid is fairly
straight-forward and can be found, for example, in Chapter 11 of McQuarrie,
2000.

   D.A. McQuarrie. *Statistical Mechanics*. University Science Books, 2000.

The treatment implemented in the :class:`CrystalThermo` class depends on
introducing normal coordinates to the entire crystal and treating each atom
in the lattice as an independent harmonic oscillator. This yields the
partition function

.. math ::
   Z = \prod_{j=1}^\text{3N} \left( \frac{e^{-\frac{1}{2}\hbar\omega_j\beta}}{1 - e^{-\hbar\omega_j\beta}} \right) e^{-E_\text{elec} \beta}

where :math:`\omega_j` are the :math:`3N` vibrational frequencies,
:math:`E_\text{elec}` is the electronic energy of the crystalline solid, and
:math:`\beta = \frac{1}{k_\text{B} T}`. Now, taking the logarithm of the
partition function and replacing the resulting sum with an integral (assuming
that the energy level spacing is essentially continuous) gives

.. math ::
   -\ln Z = E_\text{elec}\beta + \int_0^\infty \left[ \ln \left( 1 - e^{-\hbar\omega\beta} \right) + \frac{\hbar\omega\beta}{2} \right]\sigma (\omega) \text{d}\omega

Here :math:`\sigma (\omega)` represents the degeneracy or phonon density of
states as a function of vibrational frequency. Once this function has been
determined (i.e. using the :mod:`ase.phonons` module), it is a simple matter to
calculate the canonical ensemble thermodynamic quantities; namely the
internal energy, the entropy and the Helmholtz free energy.

.. math ::
   U(T) &= -\left( \frac{\partial \ln Z}{\partial \beta} \right)_\text{N,V} \\
        &= E_\text{elec} + \int_0^\infty \left[ \frac{\hbar \omega}{e^{\hbar \omega \beta} - 1} + \frac{\hbar \omega}{2} \right]\sigma (\omega) \text{d}\omega

.. math ::
   S(T) &= \frac{U}{T} + k_\text{B} \ln Z \\
        &= \int_0^\infty \left[ \frac{\hbar \omega}{T} \frac{1}{e^{\hbar \omega \beta} - 1} - k_\text{B} \ln \left(1 - e^{-\hbar \omega \beta} \right) \right]\sigma (\omega) \text{d}\omega

.. math ::
   F(T) = U(T) - T\, S(T,P)
