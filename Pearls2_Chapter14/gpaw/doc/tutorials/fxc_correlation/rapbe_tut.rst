.. _rapbe_tut:

===============================
Correlation energies from TDDFT
===============================

The Random Phase Approximation (RPA) for correlation energies comprises a nice non-empirical expression 
for the correlation energy that can be naturally combined with exact exchange to calculate binding energies. 
Due to the non-local nature of the approximation, RPA gives a good account of dispersive forces and is the 
only xc approximation capable of describing intricate binding regimes where covalent and van der Waals 
interactions are equally important. However, RPA does not provide a very accurate description of strong 
covalent bonds and typically performs worse than standard GGA functionals for molecular atomization 
energies and cohesive energies of solids.

In general, the exact correlation energy can be expressed in terms of the exact response function as 

.. math::

  E_c = -\int_0^{\infty}\frac{d\omega}{2\pi}\int_0^1d\lambda\text{Tr}\Big[v\chi^\lambda(\omega)-v\chi^{KS}(\omega)\Big],

and the RPA approximation for correlation energies is simply obtained from the RPA approximation for the response function. 
From the point of view of TDDFT, the response function can be expressed exactly in terms of the Kohn-Sham response 
function and the exchange-correlation kernel :math:`f_{xc}`:

.. math::

  \chi^\lambda(\omega) = \chi^KS(\omega) + \chi^KS(\omega)\Big[\lambda v+f_{xc}^\lambda(\omega)\Big]\chi^\lambda(\omega).

The RPA is obtained by neglecting the exchange-correlation kernel and it should be possible to improve RPA by including 
simple approximations for this kernel. However, if one tries to use a simple adiabatic kernel, one encounters severe 
convergence problems and results become worse than RPA! The reason for this is that the locality of adiabatic kernels 
renders the pair-correlation function divergent. As it turns out, the adiabatic correlation hole can be renormalized 
by a simple non-empirical procedure, which results in a density-dependent non-locality in the kernel. This can be done 
for any adiabatic kernel and the method has implemented for LDA and PBE. We refer to these approximations as renormalized 
adiabatic LDA (rALDA) and renormalized adiabatic PBE (rAPBE). We only include the exchange part of the kernel, since this 
part is linear in :math:`\lambda` and the kernel thus only needs to be evaluated for :math:`\lambda=1`.

For more details on the theory and implementation of RPA we refer to :ref:`rpa` and the tutorial :ref:`rpa_tut`. 
The RPA tutorial should be studied before the present tutorial, which inherits much of the terminology from RPA. 
Details on the theory, implementation and benchmarking of the renormalized kernels can be found in Refs. [#Olsen1]_, [#Olsen2]_, and [#Olsen3]_.

Below we give examples on how to calculate the correlation energy of a Hydrogen atom as well as the rAPBE atomization energy 
of a :math:`CO` molecule and the rAPBE cohesive energy of diamond. 
Note that some of the calculations in this tutorial will need a lot of CPU time and are essentially not possible without a supercomputer.

Example 1: Correlation energy of the Hydrogen atom
==================================================

As a first demonstration of the deficiencies of RPA, we calculate the correlation energy of a Hydrogen atom. 
The exact correlation energy should vanish for any one-electron system. The calculations can be performed with the following scripts,
starting with a standard DFT-LDA calculation:

.. literalinclude:: H.ralda_01_lda.py

followed by an RPA calculation:

.. literalinclude:: H.ralda_02_rpa_at_lda.py

and finally one using the rALDA kernel:

.. literalinclude:: H.ralda_03_ralda.py

The analogous set of scripts for PBE/rAPBE are :download:`H.ralda_04_pbe.py`, :download:`H.ralda_05_rpa_at_pbe.py`
and :download:`H.ralda_06_rapbe.py`.
The computationally-heavy RPA/rALDA/rAPBE parts can be parallelized efficiently on multiple CPUs. 
After running the scripts the LDA and PBE correlation energies may be found in the file ``H.ralda.DFT_corr_energies.txt``.
Note that a rather small unit cell is used and the results may not be completely converged with respect 
to cutoff and unit cell. Also note that the correlation energy is calculated at different cutoff energies up to 
300 eV and the values based on two-point extrapolation is printed at the end (see :ref:`rpa_tut` and :ref:`rpa` for a 
discussion on extrapolation). The results in eV are summarized below.

=====   =======   ======
 LDA    RPA       rALDA 
=====   =======   ======
-0.56    -0.55    -0.029
=====   =======   ======

=====   =======    ======
PBE     RPA        rAPBE
=====   =======    ======
-0.16    -0.55     -0.007
=====   =======    ======

The fact that RPA gives such a dramatic underestimation of the correlation energy is a general problem with the method, 
which is also seen for bulk systems. For example, for the homogeneous electron gas RPA underestimates the correlation energy 
by ~0.5 eV per electron for a wide range of densities. 
 
Example 2: Atomization energy of CO
===================================

Although RPA severely underestimates absolute correlation energies in general, energy differences are often of decent quality 
due to extended error cancellation. Nevertheless, RPA tends to underbind and performs slightly worse than PBE for atomization 
energies of molecules. The following example shows that rAPBE not only corrects the absolute correlation energies, but also 
performs better than RPA for atomization energies.

First we set up a ground state calculation with lots of unoccupied bands. This is done with the script:

.. literalinclude:: CO.ralda_01_pbe+exx.py

which takes on the order of 6-7 CPU hours. The script generates three gpw files containing the wavefunctions,
which are the input to the rAPBE calculation. The PBE and non-selfconsistent Hartree-Fock atomization energies 
are also calculated and written to the file ``CO.ralda.PBE_HF_CO.dat``. 
Next we calculate the RPA and rAPBE energies for CO with the script

.. literalinclude:: CO.ralda_02_CO_rapbe.py

The energies for C and O are obtained from the corresponding scripts
:download:`CO.ralda_03_C_rapbe.py` and :download:`CO.ralda_04_O_rapbe.py`.
The results for various cutoffs are written to the files like ``CO.ralda_rpa_CO.dat``
and ``CO.ralda_rapbe_CO.dat``.
We also print the correlation energies of the C atom to be used in a tutorial below. 
As in the case of RPA the converged result is obtained by extrapolation using the script 

.. literalinclude:: CO.ralda_05_extrapolate.py

If pylab is installed, the plot=False can be change to plot=True to visualize the quality of the extrapolation. The final results are displayed below

======   =====   =====   ======       ============
PBE      HF      RPA     rAPBE        Experimental
======   =====   =====   ======       ============
11.71    7.36    10.60    11.31         11.23
======   =====   =====   ======       ============

Example 3: Cohesive energy of diamond
=====================================
The error cancellation in RPA works best when comparing systems with similar electronic structure. 
In the case of cohesive energies of solids where the bulk energy is compared to the energy of isolated atoms, 
RPA becomes even worse than for atomization energies of molecules. Here we illustrate this for the cohesive 
energy of diamond and show that the rAPBE approximation corrects the problem. 
The initial orbitals are obtained with the script

.. literalinclude:: diamond.ralda_01_pbe.py

which takes roughly 5 minutes. The script generates ``diamond.ralda.pbe_wfcs.gpw``
and uses a previous calculation of the C atom to calculate the EXX and PBE cohesive 
energies that are written to ``diamond.ralda.PBE_HF_diamond.dat``. The RPA and rAPBE 
correlation energies are obtained with the script:

.. literalinclude:: diamond.ralda_02_rapbe_rpa.py

This takes on the order of 30 CPU hours, but can be parallelized efficiently. Finally the 
correlation part of the cohesive energies are obtained by extrapolation with the script

.. literalinclude:: diamond.ralda_03_extrapolate.py

The results are summarized below

====   ====   ====   ======       ============
PBE     HF     RPA    rAPBE       Experimental
====   ====   ====   ======       ============
7.75   5.17   7.04     7.61             7.55
====   ====   ====   ======       ============

As anticipated, RPA severly underestimates the cohesive energy, while PBE performs much better, and rAPBE comes very close to the experimental value.

.. [#Olsen1] T. Olsen and K. S. Thygesen
              *Phys. Rev. B* **86**, 081103(R) (2012)

.. [#Olsen2] T. Olsen and K. S. Thygesen
              *Phys. Rev. B* **88**, 115131 (2013)

.. [#Olsen3] T. Olsen and K. S. Thygesen
              *Phys. Rev. Lett* **112**, 203001 (2014)
