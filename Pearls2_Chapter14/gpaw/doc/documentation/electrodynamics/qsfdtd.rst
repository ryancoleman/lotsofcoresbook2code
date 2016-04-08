.. _qsfdtd:

================================================
Quasistatic Finite-Difference Time-Domain method
================================================
The optical properties of all materials depend on how they
respond (absorb and scatter) to external electromagnetic fields.
In classical electrodynamics, this response is described by
the Maxwell equations. One widely used method for solving them
numerically is the finite-difference time-domain (FDTD)
approach. \ [#Taflove]_.
It is based on propagating the electric and magnetic fields
in time under the influence of an external perturbation (light)
in such a way that the observables are expressed in real space
grid points. The optical constants are obtained by analyzing
the resulting far-field pattern. In the microscopic limit of
classical electrodynamics the quasistatic approximation is
valid and an alternative set of time-dependent equations for
the polarization charge, polarization current, and the
electric field can be derived.\ [#coomar]_

The quasistatic formulation of FDTD is implemented in GPAW.
It can be used to model the optical properties of metallic
nanostructures (i) purely classically, or (ii) in combination with
:ref:`timepropagation`, which yields :ref:`hybridscheme`.

.. TODO: a schematic picture of classical case and hybrid case

-------------------------
Quasistatic approximation
-------------------------
The quasistatic approximation of classical electrodynamics
means that the retardation effects due to the finite speed
of light are neglected. It is valid at very small length
scales, typically below ~50 nm.

Compared to full FDTD, quasistatic formulation has some
advantageous features. The magnetic field is negligible
and only the longitudinal electric field need to be
considered, so the number of degrees of freedom is
smaller. Because the retardation effects and
propagating solutions are excluded, longer time steps
and a simpler treatment of the boundary conditions can
be used.

------------
Permittivity
------------

In the current implementation, the permittivity of the classical material is parametrized as

.. math::

    \epsilon(\mathbf{r}, \omega) = \epsilon_{\infty} + \sum_j \frac{\epsilon_0 \beta_j(\mathbf{r})}{\bar{\omega}_j^2(\mathbf{r})-\mbox{i}\omega\alpha_j(\mathbf{r})-\omega^2},

where :math:`\alpha_j, \beta_j, \bar{\omega}_j` are
fitted to reproduce the experimental permittivity.
For gold and silver they can be found in Ref. \ [#Coomar]_.
Permittivity defines how classical charge density polarizes
when it is subject to external electric fields.
The time-evolution for the charges in GPAW is performed with
the leap-frog algorithm, following Ref. \ [#Gao]_.

-------------------
Geometry components
-------------------
Several routines are available to generate the basic shapes:

* :math:`\text{PolarizableBox}(\mathbf{r}_1, \mathbf{r}_2, \epsilon({\mathbf{r}, \omega}))` where :math:`\mathbf{r}_1` and :math:`\mathbf{r}_2` are the corner points, and :math:`\epsilon({\mathbf{r}, \omega})` is the permittivity inside the structure
* :math:`\text{PolarizableSphere}(\mathbf{p}, r, \epsilon({\mathbf{r}, \omega}))` where :math:`\mathbf{p}` is the center and :math:`r` is the radius of the sphere
* :math:`\text{PolarizableEllipsoid}(\mathbf{p}, \mathbf{r}, \epsilon({\mathbf{r}, \omega}))` where :math:`\mathbf{p}` is the center and :math:`\mathbf{r}` is the array containing the three radii
* :math:`\text{PolarizableRod}(\mathbf{p}, r, \epsilon({\mathbf{r}, \omega}), c)` where :math:`\mathbf{p}` is an array of subsequent corner coordinates, :math:`r` is the radius, and :math:`c` is a boolean denoting whether the corners are rounded
* :math:`\text{PolarizableTetrahedron}(\mathbf{p}, \epsilon({\mathbf{r}, \omega}))` where :math:`\mathbf{p}` is an array containing the four corner points of the tetrahedron



These routines can generate many typical geometries, and for general cases a set of tetrahedra can be used.

----------------
Optical response
----------------
The QSFDTD method can be used to calculate the optical photoabsorption
spectrum just like in :ref:`timepropagation`:
The classical charge density is first perturbed with an instantaneous
electric field, and then the time dependence of the induced dipole moment
is recorderd. Its Fourier transformation gives the photoabsorption spectrum.

-------------------------------------------
Example: photoabsorption of gold nanosphere
-------------------------------------------
This example calculates the photoabsorption spectrum of a nanosphere
that has a diameter of 10 nm, and compares the result with analytical
Mie scattering limit.

.. literalinclude:: gold_nanosphere_calculate.py

Here the *QSFDTD* object generates a dummy quantum system that is treated using
GPAW in *qsfdtd.ground_state*. One can pass the GPAW
arguments, like *xc* or *nbands*, to this function: in the example
script one empty KS-orbital was included (*nbands* =1) because GPAW
needs to propagate something. Similarly, the arguments for TDDFT
(such as *propagator*) can be passed to *time_propagation* method.

Note that the permittivity was initialized as PermittivityPlus, where
Plus indicates that a renormalizing Lorentzian term is included; this extra
term brings the static limit to vacuum value, i.e.,
:math:`\epsilon(\omega=0)=\epsilon_0`, see Ref. \ [#Sakko]_ for
detailed explanation.

The above script generates the photoabsorption spectrum and compares
it with analytical formula of the Mie theory:

.. math::
    S(\omega) = \frac{3V\omega}{2\pi^2}\mbox{Im}\left[\frac{\epsilon(\omega)-1}{\epsilon(\omega)+2}\right],

where *V* is the nanosphere volume:

|qsfdtd_vs_mie|

.. |qsfdtd_vs_mie| image:: qsfdtd_vs_mie.png

The general shape of Mie spectrum, and especially the
localized surface plasmon resonance (LSPR) at 2.5 eV,
is clearly reproduced by QSFDTD. The shoulder
at 1.9 eV and the stronger overall intensity are examples of
the inaccuracies of the used discretization scheme: the shoulder
originates from spurious surface scattering, and the intensity
from the larger volume of the nanosphere defined in the grid.
For a better estimate of the effective volume, you can take
a look at the standard output where the "Fill ratio" tells that
18.035% of the grid points locate inside the sphere. This
means that the volume (and intensity) is roughly 16% too large:

 :math:`\frac{V}{V_{\text{sphere}}}\approx\frac{0.18035\times(15\text{nm})^3)}{\frac{4}{3}\pi\times(5\text{nm})^3}\approx1.16`.

-----------
Limitations
-----------

* The scattering from the spurious surfaces of materials, which
  are present because of the representation of the polarizable
  material in uniformly spaced grid points, can cause unphysical
  broadening of the spectrum.
* Nonlinear response (hyperpolarizability) of the classical
  material is not supported, so do not use too large external
  fields. In addition to nonlinear media, also other special
  cases (nonlocal permittivity, natural birefringence, dichroism,
  etc.) are not enabled.
* The frequency-dependent permittivity of the classical material must be
  represented as a linear combination of Lorentzian oscillators. Other
  forms, such as Drude terms, should be implemented in the future. Also,
  the high-frequency limit must be vacuum permittivity. Future
  implementations should get rid of also this limitation.
* Only the grid-mode of GPAW (not e.g. LCAO) is supported.

-----------------
Technical remarks
-----------------

* Double grid technique: the calculation always uses two grids:
  one for the classical part and one for the TDDFT part. In
  purely classical simulations, suchs as the ones discussed in
  this page, the quantum subsystem contains one empty Kohn-Sham
  orbital. For more information, see the description of
  :ref:`hybridscheme` because there the double grid is very important.
* Parallelizatility: QSFDTD calculations can by parallelized
  only over domains, so use either *communicator=serial_comm* or
  *communicator=world* when initializing *QSFDTD* (or
  *FDTDPoissonSolver*) class. The domain parallelization of
  QSFDTD does not affect the parallelization of DFT calculation.
* Multipole corrections to Poissonsolver: QSFDTD module is mainly
  intended for nanoplasmonic simulations. There the charge oscillations
  are strong and the usual zero boundary conditions for the
  electrostatic potential can give inaccurate results if the simulation
  box is not large enough. In some cases, such as for single nanospheres,
  one can improve the situation by defining remove_moments argument in
  FDTDPoissonSolver: this will then use the multipole moments correction
  scheme, see e.g. Ref. \ [#Castro]_.



----
TODO
----

* Dielectrics (:math:`\epsilon_{\infty}\neq\epsilon_0`)
* Geometries from 3D model files
* Subcell averaging
* Full FDTD (retardation effects) or interface to an external FDTD software

----------------------
Combination with TDDFT
----------------------
The QSFDTD module is mainly aimed to be used in combination with :ref:`timepropagation`:
see :ref:`hybridscheme` for more information.


----------
References
----------

.. [#Taflove] A. Taflove and S. Hagness,
            Computational Electrodynamics: The Finite-Difference Time-Domain Method (3rd ed.),
            Artech House, Norwood, MA (2005).

.. [#Coomar] A. Coomar, C. Arntsen, K. A. Lopata, S. Pistinner and D. Neuhauser,
            Near-field: a finite-difference time-dependent method for simulation of electrodynamics on small scales,
            *J. Chem. Phys.* **135**, 084121 (2011)

.. [#Gao] Y. Gao and D. Neuhauser,
            Dynamical quantum-electrodynamics embedding: Combining time-dependent density functional theory and the near-field method
            *J. Chem. Phys.* **137**, 074113 (2012)

.. [#Sakko] A. Sakko, T. P. Rossi and R. M. Nieminen,
            Dynamical coupling of plasmons and molecular excitations by hybrid quantum/classical calculations: time-domain approach
            *J. Phys.: Condens. Matter* **26**, 315013 (2014)

.. [#Castro] A. Castro, A. Rubio, and M. J. Stott
             Solution of Poisson's equation for finite systems using plane wave methods
             *Canad. J. Phys.:* **81**, 1151 (2003)
