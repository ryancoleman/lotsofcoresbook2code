.. _exercise rpa:

============================================
RPA calculation of the cohesive energy of Si
============================================

In this exercise we will use GPAW to calculate the cohesive energy of 
silicon.  When calculating total energies, we shall split the 
exchange-correlation energy into an "exact" exchange contribution and 
calculate the remaining correlation energy in the random-phase 
approximation (RPA).  A comprehesive introduction to the RPA can be 
found in [Ren]_, and you are also advised to look at :ref:`this page <rpa>`.

.. note::
    
    RPA calculations are typically rather time-consuming, so in 
    this exercise we shall cut a few corners!

We obtain the cohesive energy `E_\text{coh}` of silicon as

.. math::

   E_\text{coh} = E_\text{at} - 0.5 E_\text{bulk}

where `E_\text{at}` is the total energy of an isolated silicon atom, and 
`E_\text{bulk}` is the total energy per primitive unit cell of bulk silicon.
Do you know where the factor of 0.5 comes from?

In DFT, we partition the total energy as

.. math::

   E_\text{DFT} = E_\text{Kin} + E_\text{ie} + E_\text{ii} + E_\text{Hartree} + E_\text{XC}

The purpose of this exercise is to explore different approximations for 
the exchange-correlation energy, `E_\text{XC}`.


PBE cohesive energy - bulk
==========================

We start with a generalized-gradient approximation to `E_\text{XC}` using the 
PBE functional:

.. math::

   E_\text{XC} = E_\text{PBE}

The script :download:`si.pbe.py` calculates the total 
energy of bulk silicon. The only parameters we need to choose are the 
plane-wave cutoff (i.e. number of plane waves in our basis set to describe 
the electron wavefunctions), the k-point sampling of the Brillouin zone, 
and the lattice parameter of Si.  We could calculate the 
lattice parameter ourselves from first principles, but in order to compare 
to previous calculations, we instead choose the experimental value of 
5.421 Ang [Harl]_.

Make sure you understand what the script is doing, and then try running
it.  Note that the calculation is not very heavy (it should take less
than a minute on a single CPU).  This reflects the fact that PBE is
a (semi)-local functional of the density, and also that even with
a high plane-wave cutoff, the small unit cell means we do not end up
requiring many plane-waves to describe the system.


PBE cohesive energy - atom
==========================

To complete the calculation of the cohesive energy, we need the
total energy of an isolated Si atom.  In contrast to the bulk case,
these calculations are expensive even with the PBE functional.  The
reason is that our calculations use periodic boundary conditions, leading
to unphysical interactions between replicas of the isolated atom.  We can
effectively remove this interaction placing the atom in a cubic unit cell of
side length `L`, and increase `L` until the replicas no longer see each other.
Unfortunately, the larger the value of `L`, the more plane waves we have,
which slows down the calculation considerably.

Therefore, for the purpose of this exercise, we shall not actually perform the 
calculations on the isolated Si atom - instead just provide the numbers as 
reference data.  In the next section a sample script will be given to show how to
generate the following numbers:

========  ===================
`L(\AA)`  `E_\text{PBE}` (eV)
========  ===================
6.0       -0.664402266578
7.0       -0.778484948334
8.0       -0.82500272946
9.0       -0.841856681349
10.0      -0.848092042293
11.0      -0.850367362642
12.0      -0.85109735188
========  ===================

The first column gives the side length (in Angstroms) of the simulation cell 
containing the isolated atom, and the second gives the total
energy in eV.

From the above data and your own calculations, calculate the cohesive energy 
of silicon using the PBE functional to describe exchange-correlation.  
Compare your result to the value of 4.55 eV reported in 
[Olsen]_.

.. note::
    
    The total energy delivered by GPAW is not an absolute value, but rather
    given with respect to a
    reference energy. It turns out that in this case, the reference
    energies cancel when calculating the cohesive energy, so we can forget
    about it here.  If in doubt, you can look for the line
    "reference = ...  in the GPAW output file.


EXX\@PBE cohesive energy - bulk
===============================

We now try a different approximation for the exchange-correlation energy,
which is

.. math::
  E_\text{XC} = E_\text{EXX}

An expression for the exact exchange energy `E_\text{EXX}` can be found e.g. in 
equation (9) of [Olsen]_.  The main points to note are that:

* it is fully nonlocal - to get the energy we must integrate over `\mathbf{r}`
  and `\mathbf{r}'`, which is expensive.  

* it requires knowledge of the wavefunctions, not just
  the density, which again makes it more expensive to compute.  

* in the formalism used here we calculate `E_\text{EXX}` non-self-consistently; 
  that is, we use one approximation for the exchange-correlation energy 
  (PBE) to obtain the wavefunctions, then use these wavefunctions to 
  construct the exchange energy under a different
  approximation.  As a result, this method is described as EXX\@PBE; had we
  used LDA to obtain the wavefunctions, we would have EXX\@LDA etc.

* How might a self-consistent calculation of the exchange energy compare
  to the Hartree-Fock method?

In order to obtain `E_\text{EXX}` from GPAW, we need to import the ``EXX`` class
from ``exx.py`` in our script.  The ``calculate`` method performs the
calculation of the exchange energy, while the ``get_total_energy`` method
returns the total energy of our system with `E_\text{XC}=E_\text{EXX}`.

The script :download:`si.pbe+exx.py` calculates the total 
energy of bulk Si in the EXX\@PBE approximation.  The calculation 
proceeds in two parts - first, a PBE calculation which is identical 
to that of the previous section.  Second, the exchange
part.  This part is much slower, and it is a good idea to run on a few
processors - it takes about 5 minutes on 4 CPUs.

The output file ``si.pbe+exx.exx_output.txt`` gives the details of the exchange
calculation and a breakdown of the exchange energy in terms of the
contributions from the core and valence electrons.  However for the purpose
of calculating the cohesive energy the quantity returned by the
``get_total_energy`` method and printed in ``si.pbe+exx.results.txt`` is more useful.


EXX\@PBE cohesive energy - atom
===============================

As before, we also need the energy of the isolated atom.  Look at (but don't
run!) the script :download:`atom/si.atom.pbe+exx.py`, which returns the
following output in ``pbe_and_exx_energies.txt``::

  #Box_side_length(A) PBE_total_energy(eV) PBE+EXX_total_energy(eV)
  6.0 -0.665810338359 9.88551793188
  7.0 -0.779861449204 9.79892076652
  8.0 -0.825944184466 9.76642864072
  9.0 -0.843144851642 9.75592425952
  10.0 -0.849110419847 9.7528049568
  11.0 -0.851370368753 9.7518000647
  12.0 -0.852243293624 9.75141580104
  13.0 -0.852570610869 9.75125973558

Note that :download:`atom/si.atom.pbe+exx.py` also contains 
some additional tweaking not required for the bulk calculation, 
most importantly spin-polarization; by Hund's
rules, we expect a spin-polarized atom to be more stable than the
non-spin-polarized case.

You now have enough information to calculate the cohesive energy in
the EXX\@PBE approximation.  Compare your value to that of 2.82 eV given in
[Olsen]_.  This number is dramatically different to 
the experimental value of 4.68 eV, and highlights the danger of 
neglecting correlation in solids!


(RPA+EXX)\@PBE cohesive energy - bulk
=====================================

Finally, we calculate `E_\text{XC}` including the correlation energy in the RPA:

.. math::
  E_\text{XC} = E_\text{EXX} + E_\text{RPA}

An expression for `E_\text{RPA}` is given as equation (8) in [Olsen]_.

The main ingredient here is the response function `\chi_0`, which is nonlocal,
energy dependent and constructed from a sum of an infinite number of
unoccupied electronic states.  Therefore like GW calculations, RPA
calculations are expensive to perform.  We also note that, like for exact
exchange, we construct `\chi_0` non-self-consistently, here using the
wavefunctions and eigenvalues obtained with the PBE functional.

The good news however is that compared to exact exchange calculations,
the RPA correlation energy tends to converge faster with respect to the number
of k-points and also the number of plane waves used to describe `\chi_0`, so we
can use a lighter computational setup.
Furthermore, there exists an empirical fix to the problem of the unoccupied
states which turns out to work rather well (more details below).

Like for exact exchange, the first part of our RPA calculation is performed
at the PBE level to obtain the ground state density.  We then use this density
to obtain the wavefunctions both of the occupied and some of the unoccupied
states.  The script :download:`si.rpa_init_pbe.py` performs 
this step; note it is essentially identical to 
:download:`si.pbe.py` apart from the all-important
``diagonalize_full_hamiltonian`` line.  However note that we have reduced
the k-point grid to a 4x4x4 sampling.

Having performed this step (which should take ~1 minute on 4 CPUs) we now
calculate the correlation energy using :download:`si.rpa.py`, 
which imports the ``RPACorrelation`` class from ``rpa.py``.  All the 
computational details are read from the ``bulk.gpw`` file; the only input 
we need specify is the number of plane waves used to describe `\chi_0`.  
Here we give a list of values, which means that the correlation energy 
will be calculated for each plane-wave cutoff (in eV).  The reason for 
this procedure is described below.  Note that in principle we also need 
to specify the number of unoccupied bands used in the construction of 
`\chi_0` - however here this choice is made by the code,
and sets the number of bands to equal the number of plane waves describing `\chi_0`.
Now, run :download:`si.rpa.py` (4 minutes, 4 CPUs).

Studying the output file ``si.rpa.rpa_output.txt``, we see that the code calculates
the contribution from each q point sequentially.  In fact by specifying the
``filename`` attribute of the ``RPACorrelation`` object we can generate a
restart file which allows GPAW to pick up from an interrupted calculation.
Once the contributions from all the q points have been calculated, they are
summed together with the appropriate q-point weights to construct the
correlation energy.  The correlation energy for each plane-wave cutoff is
printed near the end of the output file, under ``Total correlation energy``.
You should see that changing the plane wave cutoff from 80 to 164 eV changes
the correlation energy by over 1 eV.


(RPA+EXX)\@PBE cohesive energy - convergence
============================================

In order to converge the correlation energy, we should increase the plane-wave
cutoff describing `\chi_0` (and implicitly, the number of empty states).
However it is noted in [Harl]_ that for the electron 
gas, one expects the correlation energy to scale as

.. math::
  E_\text{RPA}(E_{cut}) = E_\text{RPA}(\infty) + A E_{cut}^{-1.5}

where `E_{cut}` is the plane-wave cutoff describing `\chi_0`.  Empirically, this
expression seems to work beyond the electron gas. 

Test this expression for silicon by plotting the correlation energy against
`E_{cut}^{-1.5}`; the intercept of the straight line should give
`E_\text{RPA}(\infty)`.  GPAW tries to guess this intercept by extrapolating
straight lines between pairs of points, and outputs the result under
``Extrapolated energies``.  How do they compare to your result?


(RPA+EXX)\@PBE cohesive energy - atom
=====================================

The corresponding scripts for the isolated atom are
:download:`atom/si.atom.rpa_init_pbe.py` and
:download:`atom/si.atom.rpa.py`. Note how, thanks to the large simulation
cell, we end up requiring almost  10000 bands for the calculation; that's a
lot of effort for a single atom!   The reference output file is
:download:`atom/si.atom.rpa_output.txt`.  Use the  data in this output file
to obtain the extrapolated correlation energy  for the single atom.

Combining the correlation energies with the EXX\@PBE calculations of the
previous section, you should now be able to calculate the cohesive energy
of silicon using exact exchange and the RPA correlation energy.  

* Compare the result of using the extrapolated correlation energies with that
  at a fixed cutoff of 164 eV.

* Compare your value to that of 4.32 eV given in [Olsen]_ 
  and the experimental value, 4.68 eV.


Conclusions
===========

After all that work, it seems that the method that gave us the cohesive
energy closest to experiment turns out to be the simplest we tried - the
generalized-gradient PBE functional.  Indeed, according to table VII of
[Harl]_, PBE outperforms EXX and RPA for a wide
range of materials.  The strength of the RPA lies in its ability to  describe
long-range correlation effects, e.g. in systems exhibiting van der Waals bonds.
Unfortunately, the complexity of these systems does not allow us to study them
in a quick exercise like this one.  Nonetheless the procedure of calculating the 
total energy employed above is exactly the same when applied to more complicated 
systems.

In order to get a consistent, high-quality description of both long-range 
and short-range correlation it is desirable to move beyond the RPA - 
but that's another story...


References
==========

.. [Ren] Ren et al., J. Mater. Sci. 47, 7447 (2012)
.. [Harl] Harl, Schimka and Kresse, Phys. Rev. B 81, 115126 (2010)
.. [Olsen] Olsen and Thygesen, Phys. Rev. B 87, 075111 (2013)
