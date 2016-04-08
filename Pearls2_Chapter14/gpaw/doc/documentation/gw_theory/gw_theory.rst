.. _gw_theory:

=======================================================
Quasi-particle spectrum in the GW approximation: theory
=======================================================

The foundations of the GW method are described in Refs. \ [#Hedin1965]_ and \ [#Hybertsen1986]_.
The implementation in GPAW is documented in Ref. \ [#Hueser2013]_.

For examples, see see :ref:`gw_tutorial`.

Introduction
============


Quasi-particle energies are obtained by replacing the DFT exchange-correlation contributions by the GW self energy and exact 
exchange:

.. math:: E_{n \mathbf{k}} = \epsilon_{n \mathbf{k}} + Z_{n \mathbf{k}} \cdot \text{Re} \left(\Sigma_{n \mathbf{k}}^{\vphantom{\text{XC}}} + \epsilon^{\text{EXX}}_{n \mathbf{k}} - V^{\text{XC}}_{n \mathbf{k}} \right)

where :math:`n` and :math:`\mathbf{k}` are band and k-point indices, respectively.

The different contributions are:

:math:`\epsilon_{n \mathbf{k}}`: Kohn-Sham eigenvalues taken from a groundstate calculation

:math:`V^{\text{XC}}_{n \mathbf{k}}`: DFT exchange-correlation contributions extracted from a groundstate calculation

:math:`\epsilon^{\text{EXX}}_{n \mathbf{k}}`: exact exchange contributions

The renormalization factor is given by:

.. math:: Z_{n \mathbf{k}} = \left(1 - \text{Re}\left< n \mathbf{k}\middle| \frac{\partial}{\partial\omega} \Sigma(\omega)_{|\omega = \epsilon_{n \mathbf{k}}}\middle| n \mathbf{k}\right>\right)^{-1}

:math:`\left| n \mathbf{k} \right>` denotes the Kohn-Sham wavefunction which is taken from the groundstate calculation.

The self energy is expanded in plane waves, denoted by :math:`\mathbf{G}` and :math:`\mathbf{G}'`:

.. math:: \Sigma_{n \mathbf{k}} =& \left<n \mathbf{k} \middle| \Sigma(\omega) \middle|n \mathbf{k} \right>\\
 =& \frac{1}{\Omega} \sum\limits_{\mathbf{G} \mathbf{G}'} \sum\limits_{\vphantom{\mathbf{G}}\mathbf{q}}^{1. \text{BZ}} \sum\limits_{\vphantom{\mathbf{G}}m}^{\text{all}} \frac{i}{2 \pi} \int\limits_{-\infty}^\infty\!d\omega'\, W_{\mathbf{G} \mathbf{G}'}(\mathbf{q}, \omega') \, \cdot \\
 & \frac{\rho^{n \mathbf{k}}_{m \mathbf{k} - \mathbf{q}}(\mathbf{G}) \rho^{n \mathbf{k}*}_{m \mathbf{k} - \mathbf{q}}(\mathbf{G}')}{\omega + \omega' - \epsilon_{m \, \mathbf{k} - \mathbf{q}} + i \eta \, \text{sgn}(\epsilon_{m \, \mathbf{k} - \mathbf{q}} - \mu)}_{|\omega = \epsilon_{n \mathbf{k}}}

where :math:`m` runs both over occupied and unoccupied bands and :math:`\mathbf{q}` covers the differences between all k-points in the first Brillouin zone. :math:`\Omega = \Omega_\text{cell} \cdot N_\mathbf{k}` is the volume and :math:`\eta` an (artificial) broadening parameter. :math:`\mu` is the chemical potential.

The screened potential is calculated from the (time-ordered) dielectric matrix in the Random Phase Approximation:

.. math:: W_{\mathbf{G} \mathbf{G}'}(\mathbf{q}, \omega) = \frac{4 \pi}{|\mathbf{q} + \mathbf{G}|} \left( (\varepsilon^{\text{RPA}-1}_{\mathbf{G} \mathbf{G}'}(\mathbf{q}, \omega) - \delta^{\vphantom{\text{RPA}}}_{\mathbf{G} \mathbf{G}'} \right) \frac{1}{|\mathbf{q} + \mathbf{G}'|}

Refer to :ref:`df_theory` for details on how the response function and the pair density matrix elements :math:`\rho^{n \mathbf{k}}_{m \mathbf{k} - \mathbf{q}}(\mathbf{G}) \equiv \left<n \mathbf{k} \middle| e^{i(\mathbf{q} + \mathbf{G})\mathbf{r}} \middle|m \, \mathbf{k} \!-\! \mathbf{q} \right>` including the PAW corrections are calculated.

Coulomb divergence
==================


The head of the screened potential (:math:`\mathbf{G} = \mathbf{G}' = 0`) diverges as :math:`1/q^2` for :math:`\mathbf{q} \rightarrow 0`. This divergence, however, is removed for an infinitesimally fine k-point sampling, as :math:`\sum\limits_{\mathbf{q}} \rightarrow \frac{\Omega}{(2\pi)^3} \int\!d^3 \mathbf{q} \propto q^2`. Therefore, the :math:`\mathbf{q} = 0` term can be evaluated analytically, which yields:

.. math:: W_{\mathbf{00}}(\mathbf{q}=0, \omega) = \frac{2\Omega}{\pi} \left(\frac{6\pi^2}{\Omega}\right)^{1/3} \varepsilon^{-1}_{\mathbf{00}}(\mathbf{q} \rightarrow 0, \omega)

for the head and similarly

.. math:: W_{\mathbf{G0}}(\mathbf{q}=0, \omega) = \frac{1}{|\mathbf{G}|} \frac{\Omega}{\pi} \left(\frac{6\pi^2}{\Omega}\right)^{2/3} \varepsilon^{-1}_{\mathbf{G0}}(\mathbf{q} \rightarrow 0, \omega)

for the wings of the screened potential. Here, the dielectric function is used in the optical limit.

This is only relevant for the terms with :math:`n = m`, as otherwise the pair density matrix elements vanish: :math:`\rho^{n \mathbf{k}}_{m \mathbf{k}} = 0` for :math:`n \neq m`.

Frequency integration
=====================
 :math:`\rightarrow` ``w = (wlin, wmax, dw)``


The frequency integration is performed numerically on a user-defined grid for positive values only. This is done by rewriting the integral as:

.. math:: & \int\limits_{-\infty}^\infty\!d\omega'\, \frac{W(\omega')}{\omega + \omega' - \epsilon_{m \, \mathbf{k} - \mathbf{q}} \pm i \eta}\\
 =& \int\limits_{0}^\infty\!d\omega'\, W(\omega') \left(\frac{1}{\omega + \omega' - \epsilon_{m \, \mathbf{k} - \mathbf{q}} \pm i \eta} + \frac{1}{\omega - \omega' - \epsilon_{m \, \mathbf{k} - \mathbf{q}} \pm i \eta}\right)

with the use of :math:`W(\omega') = W(-\omega')`.

The frequency grid is defined by an array of three values: :math:`[\omega_{\text{lin}}, \omega_{\text{max}}, \Delta\omega]`. This creates a linear grid from :math:`0` to :math:`\omega_{\text{lin}}` with a spacing :math:`\Delta\omega`. Above :math:`\omega_{\text{lin}}`, the spacing increases linearly with every step up to the maximum frequency :math:`\omega_{\text{max}}`. All values are in eV. The maximum frequency has to be bigger than the largest transition energy :math:`|\epsilon_{n \, \mathbf{k}} - \epsilon_{m \, \mathbf{k} - \mathbf{q}}|` included in the calculation. 

Plasmon Pole Approximation
==========================
 :math:`\rightarrow` ``ppa = True``


Within the plasmon pole approximation (PPA), the dielectric function is modelled as a single peak at the main plasmon frequency :math:`\tilde{\omega}_{\mathbf{G}\mathbf{G}'}(\mathbf{q})`:

.. math:: \varepsilon^{-1}_{\mathbf{G}\mathbf{G}'}(\mathbf{q}, \omega) = R _{\mathbf{G}\mathbf{G}'}(\mathbf{q}) \left(\frac{1}{\omega - \tilde{\omega}_{\mathbf{G}\mathbf{G}'}(\mathbf{q}) + i\eta} - \frac{1}{\omega + \tilde{\omega}_{\mathbf{G}\mathbf{G}'}(\mathbf{q}) - i\eta}\right)

The two parameters are found by fitting this expression to the full dielectric function for the values :math:`\omega = 0` and :math:`\omega = i E_0`:

.. math:: \varepsilon^{-1}_{\mathbf{G}\mathbf{G}'}(\mathbf{q}, 0) =& \frac{-2 R}{\tilde{\omega}} \hspace{0.5cm} \varepsilon^{-1}_{\mathbf{G}\mathbf{G}'}(\mathbf{q}, iE_0) = \frac{-2 R \tilde{\omega}}{E_0^2 + \tilde{\omega}^2}\\
 \Rightarrow \tilde{\omega}_{\mathbf{G}\mathbf{G}'}(\mathbf{q}) =& E_0 \sqrt{\frac{\varepsilon^{-1}_{\mathbf{G}\mathbf{G}'}(\mathbf{q}, iE_0)} {\varepsilon^{-1}_{\mathbf{G}\mathbf{G}'}(\mathbf{q}, 0) - \varepsilon^{-1}_{\mathbf{G}\mathbf{G}'}(\mathbf{q}, iE_0)}}\\
 R _{\mathbf{G}\mathbf{G}'}(\mathbf{q}) =& -\frac {\tilde{\omega}_{\mathbf{G}\mathbf{G}'}(\mathbf{q})}{2} \varepsilon^{-1}_{\mathbf{G}\mathbf{G}'}(\mathbf{q}, 0)

In this way, the frequency integration for the self energy can be evaluated analytically. The fitting value :math:`E_0` has to be chosen carefully. By default, it is 1 H.

Static COHSEX
==========================
 :math:`\rightarrow` ``w = None``


In the static limit :math:`\omega - \epsilon_{m \, \mathbf{k} - \mathbf{q}} = 0`, the self energy can be split into two terms, which can be identified as screened exchange and Coulomb hole:

.. math:: \Sigma_{n \mathbf{k}}^{\text{SEX}} = - \frac{1}{\Omega} \sum\limits_{\mathbf{G} \mathbf{G}'} \sum\limits_{\vphantom{\mathbf{G}}\mathbf{q}} \sum\limits_{\vphantom{\mathbf{G}}m}^{\text{occ}} \varepsilon^{-1}_{\mathbf{G} \mathbf{G}'}(\mathbf{q}, 0) V_{\mathbf{G} \mathbf{G}'}^{\vphantom{-1}}(\mathbf{q}) \rho^{n \mathbf{k}}_{m \mathbf{k} - \mathbf{q}}(\mathbf{G}) \rho^{n \mathbf{k}*}_{m \mathbf{k} - \mathbf{q}}(\mathbf{G}')

.. math:: \Sigma_{n \mathbf{k}}^{\text{COH}} = \frac{1}{2 \Omega} \sum\limits_{\mathbf{G} \mathbf{G}'} \sum\limits_{\vphantom{\mathbf{G}}\mathbf{q}} \sum\limits_{\vphantom{\mathbf{G}}m}^{\text{all}} \left(\varepsilon^{-1}_{\mathbf{G} \mathbf{G}'}(\mathbf{q}, 0) - \delta_{\mathbf{G} \mathbf{G}'}^{\vphantom{-1}}\right) V_{\mathbf{G} \mathbf{G}'}^{\vphantom{-1}}(\mathbf{q}) \rho^{n \mathbf{k}}_{m \mathbf{k} - \mathbf{q}}(\mathbf{G}) \rho^{n \mathbf{k}*}_{m \mathbf{k} - \mathbf{q}}(\mathbf{G}')

where :math:`V_{\mathbf{G} \mathbf{G}'}(\mathbf{q}) = 4 \pi / |\mathbf{q} + \mathbf{G}||\mathbf{q} + \mathbf{G}'|` is the Coulomb potential.

The quasi-particle energies are then calculated as:

.. math::  E_{n \mathbf{k}} = \epsilon_{n \mathbf{k}} + \Sigma_{n \mathbf{k}}^{\text{SEX}} + \Sigma_{n \mathbf{k}}^{\text{COH}} - V^{\text{XC}}_{n \mathbf{k}}
 
Hilbert transform
=================


Currently, there are two different methods implemented for evaluating the self energy.

Method 1 (which is the default ``hilbert_trans = False``) performs the summation over plane waves first:

.. math:: \sum\limits_{\mathbf{G} \mathbf{G}'} W_{\mathbf{G} \mathbf{G}'}(\mathbf{q}, \omega') \rho^{n \mathbf{k}}_{m \mathbf{k} - \mathbf{q}}(\mathbf{G}) \rho^{n \mathbf{k}*}_{m \mathbf{k} - \mathbf{q}}(\mathbf{G}')

Then, the frequency integration with

.. math:: \frac{1}{\omega + \omega' - \epsilon_{m \, \mathbf{k} - \mathbf{q}} + i \eta \, \text{sgn}(\epsilon_{m \, \mathbf{k} - \mathbf{q}} - \mu)} \hspace{0.4cm} \textsf{and} \hspace{0.4cm} - \frac{1}{\left(\omega + \omega' - \epsilon_{m \, \mathbf{k} - \mathbf{q}} + i \eta \, \text{sgn}(\epsilon_{m \, \mathbf{k} - \mathbf{q}} - \mu)\right)^2}

for the self energy and its derivative is carried out, where :math:`\omega = \epsilon_{n \mathbf{k}}`. This is done for every :math:`(n \, \mathbf{k})` and :math:`(m \, \mathbf{k}\!-\!\mathbf{q})` seperately.

Method 2 (``hilbert_trans = True``) reverses this order by doing the frequency integration first for all :math:`\omega` on the grid. Then, for every :math:`(n \, \mathbf{k})` and :math:`(m \, \mathbf{k}\!-\!\mathbf{q})`, the contributions to :math:`\Sigma(\omega = \epsilon_{n \mathbf{k}})` and its derivative are found by linear interpolation using the two closest points on the frequency grid with :math:`\omega_i \leq \omega = \epsilon_{n \mathbf{k}} < \omega_{i+1}`. For :math:`\omega = 0`, three points are used for the interpolation. This is similar to using the Hilbert transform for the dielectric response function.

While the first method is more accurate, the second method can reduce the computational costs significantly. As long as the chosen frequency grid is fine enough, both methods yield the same results.

See ref. \ [#Kresse2006]_ for details.

Parallelization
===============
 :math:`\rightarrow` ``wpar = int``


By default, the calculation is fully parallelized over k-points, that means all :math:`\mathbf{q}` in the summation. When more memory is required for storing the dielectric matrix as a function of :math:`\omega`, :math:`\mathbf{G}` and :math:`\mathbf{G}'`, additional parallelization over frequencies may be necessary. This can be done by increasing ``wpar``. This value determines over how many CPUs the dielectric function (and its related quantities) should be distributed. Information about the memory usage is printed in the output file ``df.out``.

Note, that ``wpar`` needs to be an integer divisor of the number of requested CPUs.

I/O
===


All necessary informations of the system are read from ``file = 'filename.gpw'`` which must contain the wavefunctions. This is done by performing ``calc.write('groundstate.gpw', 'all')`` after the groundstate calculation. GW supports grid mode and planewave basis.

Especially for big systems, it might be reasonable to determine the exact exchange contributions seperately and store them in a pickle file which can be read by defining ``exxfile = 'filename.pckl'`` (see below). The band and k-point indices must match the ones used for the GW calculation. The pickle file needs to contain the following data:

================= ==============================================================================
``gwkpt_k``       list of k-point indices

``gwbands_n``     list of bands indices

``e_skn``         DFT eigenvalues as array with spin, k-points and bands

``vxc_skn``       DFT exchange-correlation contributions as array with spin, k-points and bands

``exx_skn``       exact exchange contributions as array with spin, k-points and bands
================= ==============================================================================

See the GW tutorial for an example: :ref:`gw_tutorial`

The output is written to ``txt = 'filename.out'`` which summarizes the input and results and gives an estimation of the timing while the calculation is running. An additional file ``df.out`` is created for the calculation of the dielectric matrix.

All results are also stored in a pickle file called ``GW.pckl`` by default, which contains all data listed in the table above and addionally ``Sigma_skn``, ``Z_skn`` and ``QP_skn`` for the self energy contributions, renormalization factors and the quasi-particle bandstructure, respectively.

Convergence
===========


The results must be converged with respect to:

- the number of k-points from the groundstate calculation
    A much finer k-point sampling might be required for converging the GW results than for the DFT bandstructure.

- the number of bands included in the calculation of the self energy ``nbands``

- the planewave energy cutoff ``ecut``
    ``ecut`` and ``nbands`` do not converge independently. As a rough estimation, ``ecut`` should be around the energy of the highest included band.

- the fineness of the frequency grid ``wlin, dw``
    The grid needs to resolve the features of the DFT spectrum.

- the broadening ``eta``
    This parameter is only used for the response function and in the plasmon pole approximation. Otherwise, it is automatically set to :math:`\eta = 4 \Delta\omega`.


Parameters
==========


=================  =================  ===================  ====================================================
keyword            type               default value        description
=================  =================  ===================  ====================================================
``file``           ``str``            None                 gpw filename
                                                           groundstate calculation including wavefunctions
``nbands``         ``int``            equal to number of   Number of bands
                                      plane waves
``bands``          ``numpy.ndarray``  equal to nbands      Band indices for QP spectrum
``kpoints``        ``numpy.ndarray``  all irreducible      K-point indices for QP spectrum
                                      k-points
``e_skn``          ``numpy.ndarray``  None                 User-defined starting point eigenvalues
``eshift``         ``float``          None                 Manual shift of unoccupied bands (in eV)
``w``              ``numpy.ndarray``  None                 [wlin, wmax, dw] for defining frequency grid (in eV)
``ecut``           ``float``          150 (eV)             Planewave energy cutoff.
``eta``            ``float``          0.1 (eV)             Broadening parameter.
``ppa``            ``bool``           False                Use Plasmon Pole Approximation
``E0``             ``float``          27.2114 (eV)         Frequency for fitting PPA
``hilbert_trans``  ``bool``           False                False = method 1, True = method 2
``wpar``           ``int``            1                    Parallelization in energy
``vcut``           ``str``            None                 Coulomb truncation (currently, only '2D' supported)
``txt``            ``str``            None                 Output filename
=================  =================  ===================  ====================================================

Functions
=========

``get_exact_exchange(ecut=None, communicator=world, file='EXX.pckl')``

calculates exact exchange and Kohn-Sham exchange-correlation contributions for given ``bands`` and ``kpoints``
and stores everything in a pickle file.

In planewave mode ``ecut`` is taken from the groundstate calculation.
Otherwise, it can be chosen independently from the actual GW calculation. The value needs to be given in eV.
Note that the exact exchange and GW may converge differently with respect to ``ecut``.

``get_QP_spectrum(exxfile='EXX.pckl', file='GW.pckl')``

calculates GW quasi-particle spectrum reading exact exchange and exchange-correlation contribution from ``exxfile``
and stores all results in a pickle file.


References
==========


.. [#Hedin1965] L. Hedin,
                "New Method for Calculating the One-Particle Green's Function with Application to the Electron-Gas Problem",
                *Phys. Rev.* **139**, A796 (1965).

.. [#Hybertsen1986] M.S. Hybertsen and S.G. Louie,
                    "Electron correlation in semiconductors and insulators: Band gaps and quasiparticle energies",
                    *Phys. Rev. B* **34**, 5390 (1986).

.. [#Hueser2013] F. Hüser, T. Olsen, and K. S. Thygesen,
                 "Quasiparticle GW calculations for solids, molecules, and two-dimensional materials",
                 *Phys. Rev. B* **87**, 235132 (2013).

.. [#Kresse2006] M. Shishkin and G. Kresse,
                 "Implementation and performance of the frequency-dependent GW method within the PAW framework",
                 *Phys. Rev. B* **74**, 035101 (2006).
