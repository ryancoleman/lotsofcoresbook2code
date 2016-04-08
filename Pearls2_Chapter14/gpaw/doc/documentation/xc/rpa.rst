.. _rpa:

=======================
RPA correlation energy
=======================

The correlation energy within the Random Phase Approximation (RPA) can be written

.. math::

  E_c^{RPA} = \int_0^{\infty}\frac{d\omega}{2\pi}\text{Tr}\Big[\text{ln}\{1-\chi^0(i\omega)v\}+\chi^0(i\omega)v\Big],
 
where `\chi^0(i\omega)` is the non-interacting (Kohn-Sham) response function evaluated at complex frequencies, `\text{Tr}` is the Trace and `\it{v}` is the Coulomb interaction. The response function and Coulomb interaction are evaluated in a plane wave basis as described in :ref:`df_tutorial` and :ref:`df_theory` and for periodic systems the Trace therefore involves a summation over `\mathbf{q}`-points, which are determined from the Brillouin zone sampling used when calculating `\chi^0(i\omega)`.

The RPA correlation energy is obtained by::
    
    from gpaw.xc.rpa import RPACorrelation
    rpa = RPACorrelation(calc, txt='rpa_correlation.txt')
    E_rpa = rpa.calculate(ecut=400)

where calc is either a calculator object containing converged wavefunctions from a ground state calculation or a string reference to a .gpw file containing wavefunctions. If calc is a calculator object it should be loaded in serial since the RPA parallellization scheme is rather different from that of standard DFT calculatons. txt denotes the output file. The RPACorrelation also takes a number of optional keywords described below. The calculate() function performs the actual calculation at the cutoff energy specified by ecut (in eV). In addition the rpa calculator will calculate the correlation energy at four values for the cutoff energies up to the specified cutoff, but one can also give a list of cutoff values instead. By default, the response function is calculated with the same number of bands as the number of plane waves, but one can also specify that it should use N bands with nbands=N in the calculate() function.


Parameters
==========

=================== ================== =================== ==================================================================
keyword             type               default value       description
=================== ================== =================== ==================================================================
``nfrequencies``    ``int``            16                  Number of Gauss-legendre points used in the
                                                           integration.
``frequency_cut``   ``float``          800. (eV)           The maximum frequency is the largest frequency
                                                           included in the Gauss-Legendre integration. The integral is
                                                           always an approximation to the infinite integral, but the
                                                           max frequency determines the distribution of frequencies.
``frequency_scale`` ``float``          2.0 (eV)            The frequency scale sets the density of frequency
                                                           points near :math:`\omega = 0`.
``frequencies``     ``numpy.ndarray``  None                Specifies frequency points used to integrate the
                                                           correlation integrand.
                                                           Ex: numpy.linspace(0,20,201). If None, the Gauss-legendre
                                                           method is used.
``weights``         ``numpy.ndarray``  None                Should be used in conjunction with frequencies (e.i.
                                                           when not using the Gauss-Legendre integration). For example
                                                           np.array([0.5,1,1,...,1,1,0.5] gives a trapezoid integration
``skip_gamma``      ``bool``           False               For metals the :math:`\mathbf{q} = 0` point can give rise
                                                           to divergent contributions and it may be faster to converge
                                                           the k-point sampling if this point is excluded.
``nblocks``         ``int``            1                   **G**-vector parallelization. Default parallelization scheme is over
                                                           kpoints, spin and bands. If memory becomes an issue it can be an
                                                           advantage to use **G**-vector parallelization also.
``filename``        ``str``            None                Restart file. If calculations with k-point sampling, the
                                                           contributions from different q-points are calculated
                                                           sequentially and written to filename such that these do not have
                                                           to be recalculated when a calculation is restarted.
=================== ================== =================== ==================================================================

In addition to the usual kpoint and plane wave cutoff, the RPA correlation energy needs to be converged with respect to a plane wave cutoff in the response function (set by ecut) and the frequency integration. As it turns out, the integrand is usually  rather smooth and one can perform the integration with 8-16 (special!) Gauss-Legendre frequency points, but see the tutorial :ref:`rpa_tut` for an example of converging the frequency integration.
        
Convergence
===========

A major complication with the RPA correlation energy is that it converges very slowly with the number of unoccupied bands included in the evaluation of `\chi^0(i\omega)`. However, as described in Ref. \ [#Harl1]_ the high energy part of the response function resembles the Lindhard function, which for high energies gives a correlation energy converging as

.. math::

  E_c^{Lindhard}(E^{\chi}_{cut}) = E_c^{\infty}+\frac{A}{(E^{\chi}_{cut})^{3/2}},

where `E^{\chi}_{cut}` is cutoff energy used in the evaluation of `\chi^0`. With an external potential, the number of unoccupied bands is an additional convergence parameter, but for reproducing the scaling of the Lindhard function, it is natural to set the total number of bands equal to the number of plane waves used. Thus, to obtain a converged RPA correlation energy one should proceed in three steps.

* Perform a ground state calculation with a lot of converged unoccupied bands.
  
* Define a list of cutoff energies - typically something like [200, 225, 250, 275, 300] (eV). For each cutoff energy perform an RPA correlation energy calculation with the number bands `n` set equal to the number of plane waves defined by that cutoff energy.

* Fit the list of obtained correlation energies to `E_c^{RPA}(E) = E_c^{\infty}+A/E^{3/2}` to obtain `E_c^{\infty}=E_c^{RPA}`.

Per default, the rpa module defines a list of five cutoff energies up to the specified value and performs the extrapolation at the end of the calculation. If one is not interested in the total correlation energy, but only energy differences between similar systems, it is sometimes possible to avoid the extrapolation procedure and the rpa correlation energy can be obtained at a single point by specifying a list with one element (for example ecut=[400]).

.. [#Harl1] J. Harl and G. Kresse,
            *Phys. Rev. B* **77**, 045136 (2008)

.. [#Harl2] J. Harl and L. Schimka and G. Kresse,
            *Phys. Rev. B* **81**, 115126 (2010)
