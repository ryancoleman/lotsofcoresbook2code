.. _dscf:

===========================
Delta Self-Consistent Field
===========================

--------------------------------------------
Linear expansion Delta Self-Consistent Field
--------------------------------------------

The method of linear expansion Delta Self-Consistent Field \ [#delscf]_ 
adds the density of a specified orbital `\varphi_a(r)` to the 
total density in each step of the self-consistency cycle. The extra charge 
is usually taken from the fermi level to keep the system neutral:

.. math::

  n(r) = \sum_nf_{N-1}(T,\varepsilon_n)|\varphi_n(r)|^2 + |\varphi_a(r)|^2.

with `N` being the total number of electrons and `f_{N-1}(T,\varepsilon_n)`
is the Fermi-Dirac distribution of the `N-1` electron system . To get the 
band energy right `\varphi_a(r)` needs to be expanded in Kohn-Sham orbitals:

.. math::

  |\varphi_a\rangle = \sum_nc_{na}|\varphi_n\rangle, 
  \qquad c_{na} = \langle\varphi_n|\varphi_a\rangle

and the band energy of the orbital becomes

.. math::

  \varepsilon_a = \sum_n|c_{na}|^2\varepsilon_n.

The method is a generalization of traditional Delta Self-Consistent Field
where only the occupation numbers are modified and it will reduce to that, 
if only one (normalized) term is included in the expansion of `\varphi_a(r)`.

----------------
Simple molecules
----------------

The example below calculates the excitation energy of the 
`5\sigma\rightarrow2\pi` transition in CO. We only specify that the 
`2\pi` orbital should be occupied ([[1.0, lumo, 1]] means 1.0 electrons 
in lumo with spin 1) and the method will take the electron from highest 
occupied orbital which in this case is `5\sigma`.

The lumo is an instance of the class AEOrbital which calculates the 
expansion of the saved `2\pi` state in each iteration step.
In order to obtain the all-electron overlaps `\langle\varphi_n|2\pi\rangle` 
we need to supply the projector overlaps in addition to the 
pseudowavefunction.

Exciting the LUMO in CO::

    from ase.structure import molecule
    from gpaw import GPAW
    from gpaw import dscf

    # Ground state calculation
    #------------------------------------------------------------------

    calc = GPAW(nbands=8, h=0.2, xc='PBE', spinpol=True,
                convergence={'energy': 100,
                             'density': 100,
                             'eigenstates': 1.0e-9,
                             'bands': -1})

    CO = molecule('CO')
    CO.center(vacuum=3)
    CO.set_calculator(calc)

    E_gs = CO.get_potential_energy()

    # Obtain the pseudowavefunctions and projector overlaps of the
    # state which is to be occupied. n=5,6 is the 2pix and 2piy orbitals
    n = 5
    molecule = [0, 1]
    wf_u = [kpt.psit_nG[n] for kpt in calc.wfs.kpt_u]
    p_uai = [dict([(molecule[a], P_ni[n]) for a, P_ni in kpt.P_ani.items()])
             for kpt in calc.wfs.kpt_u]

    # Excited state calculation
    #--------------------------------------------

    calc_es = GPAW(nbands=8, h=0.2, xc='PBE', spinpol=True,
                   convergence={'energy': 100,
                                'density': 100,
                                'eigenstates': 1.0e-9,
                                'bands': -1})

    CO.set_calculator(calc_es)
    lumo = dscf.AEOrbital(calc_es, wf_u, p_uai)
    #lumo = dscf.MolecularOrbital(calc, weights={0: [0, 0, 0,  1],
    #                                            1: [0, 0, 0, -1]})
    dscf.dscf_calculation(calc_es, [[1.0, lumo, 1]], CO)

    E_es = CO.get_potential_energy()

    print 'Excitation energy: ', E_es-E_gs

The commented line ``lumo = dscf.Molecular...`` 
uses another class to specify the `2\pi` orbital of CO which does not require 
a ground state calculation of the molecule. In the simple example above the 
two methods give identical results, but for more complicated systems the 
AEOrbital class should be used \ [#des]_. When using the AEOrbital class 
a new calculator object must be constructed for the dscf calculation.

In the example above we only specify a single state, but the function 
``dscf.dscf_calculation`` takes a list of orbitals as input and we could for 
example have given the argument [[1.0, lumo, 1], [-1.0, pi, 0]] which would 
force the electron to be taken from the `\pi` orbital with spin 0. The pi 
should of course be another instance of the AEOrbital class.

---------------------
Exciting an adsorbate
---------------------
The method of linear expansion Delta Self-Consistent Field was designed
for calculations with strongly hybridized orbitals. For example molecules 
chemisorbed on transition metals. In such cases the 
traditional Delta Self-Consistent Field breaks down since the orbital
to be occupied is no longer well described by a single Kohn-Sham state.

The script :svn:`~doc/documentation/dscf/homo.py` calculates 
the HOMO energy of CO adsorbed on-top Pt(111). The script starts
from scratch, but usually one would start from an optimized configuration
saved in a file ``gs.gpw``. The script only calculates the total energy of 
the excited state so the excitation energy is obtained as the difference 
between ground and excited state energies.

First a calculation of gas-phase CO is performed and the 
HOMO pseudo-wavefunctions and the projector overlaps are saved. The 
energy range [-100.0, 0.0] means we only include states below the Fermi
level (default is states above).

The script :svn:`~doc/documentation/dscf/lumo.py` calculates
the LUMO energy of the same system, but is slightly more complicated due to 
the degeneracy of the `2\pi` orbital. We would like to occupy the `2\pi_y` 
orbital and  we need to figure out which band (5 or 6) this orbital 
corresponds to in each k-point before we start the slab calculation.

.. [#delscf] J. Gavnholt, T. Olsen, M. Engelund and J. Schiøtz,
             Delta Self-Consistent Field as a method to obtain potential
	     energy surfaces of excited molecules on surfaces,
             *Phys. Rev. B* **78**, 075441 (2008)

.. [#des]    T. Olsen, J. Gavnholt and J. Schiøtz,
             Hot electron mediated desorption rates calculated from excited
	     state potential energy surfaces,
             *Phys. Rev. B* **79**, 035403 (2009)
