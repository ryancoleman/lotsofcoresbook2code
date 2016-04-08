.. _gw_tutorial:

=========================================================
Quasi-particle spectrum in the GW approximation: tutorial
=========================================================

For a brief introduction to the GW theory and the details of its implementation in GPAW, see :ref:`gw_theory`.

More information can be found here:

    \F. Hüser, T. Olsen, and K. S. Thygesen

    `Quasiparticle GW calculations for solids, molecules, and two-dimensional materials`__

    Physical Review B, Vol. **87**, 235132 (2013)

    __ http://prb.aps.org/abstract/PRB/v87/i23/e235132


Quasi-particle spectrum of bulk silicon
=======================================


groundstate calculation
-----------------------

First, we need to do a regular groundstate calculation.
We do this in plane wave mode and choose the LDA exchange-correlation functional.
In order to keep the computational efforts small, we start with (3x3x3) k-points and a plane wave basis up to 200 eV.

.. literalinclude:: Si_groundstate.py

It takes a few minutes on a single CPU.
The last line in the script creates a .gpw file which contains all the informations of the system, including the wavefunctions.

.. note::

    You can change the number of bands to be written out by using ``calc.diagonalize_full_hamiltonian(nbands=...)``.
    This will be useful for higher plane wave cutoffs.

.. note::

    For calculations including only the :math:`\Gamma` point, insert ``dtype=complex`` in the calculator.
    This holds both for the plane wave and the grid mode.

the GW calculator
-----------------

Next, we set up the GW calculator, where we define all the required parameters
as well as the k-point and band indices for which we want to calculate the quasi-particle spectrum.
Here, we do this for the complete irreducible Brioullin zone and 4 bands around the Fermi level
(silicon has 8 valence electrons and the bands are double occupied, starting from band index 0,
so the corresponding band indices are 2,3,4 and 5).

.. literalinclude:: Si_gw.py
    :lines: 1-12

calculating the exact exchange contributions
--------------------------------------------

It is highly recommended (though not necessary) to start with calculating the exact exchange contributions.
This is simply done by calling:

.. literalinclude:: Si_gw.py
    :lines: 14

In the output file, we find the results for non-selfconsistent Hartree-Fock,
sorted by spin, k-points (rows) and bands (columns).

.. note::

    By default, the results are stored in a pickle file called ``EXX.pckl``.
    The name can be changed by using ``gw.get_exact_exchange(file='myown_EXX.pckl')``.

Check for convergence with respect to the plane wave cutoff energy and number of k points
by changing the respective values in the groundstate calculation and restarting.

.. image:: Si_EXX.png
       :height: 400 px

A k-point sampling of (9x9x9) and 150 eV plane wave cutoff seems to give well converged results.

calculating the self-energy
---------------------------

Now, we are ready to calculate the GW quasiparticle spectrum by calling:

.. literalinclude:: Si_gw.py
    :lines: 16

While the calculation is running, timing information is printed in the output file.

In the end, the results for the quasiparticle spectrum are printed,
again sorted by spin, k-points (rows) and bands (columns).

.. note::

    By default, the results are stored in a pickle file called ``GW.pckl``.
    The name can be changed by using ``gw.get_QP_spectrum(file='myown_GW.pckl', exxfile='myown_EXX.pckl')``.

Again, results need to be checked for convergence with respect to k-point sampling, plane wave cutoff energy and number of bands.

A reasonable choice are the values determined previously for the Hartree Fock band gap, that means (9x9x9) k points and plane waves up to 150 eV.
Usually, the self-energy converges faster with respect to the number of k points than the exchange contributions, whereas the dependence on the plane wave cutoff energy is similar.
The number of bands should be chosen in a way, so that the energy of the highest band is around `E_\text{cut}`. This is easiest done by using the default value ``nbands = None``,
which sets the number of bands equal to the number of plane waves.

frequency dependence
--------------------

Next, we should check the quality of the Plasmon Pole Approximation and use the fully frequency-dependent dielectric matrix.

Remove the line ``ppa = True`` and insert ``w = np.array([50., 150., 0.05])``.
This creates a frequency grid which is linear from 0 to 50 eV with a spacing of 0.05 eV and increasing steps up to 150 eV.
This will correspond to ~1064 frequency points. The calculation takes about 1 hour.

The results should be very close to what we obtained within the Plasmon Pole Approximation, verifying its validity for this system.

At last, see how the results depend on the chosen frequency grid. It is important to have a very fine grid in the lower frequency range,
where the electronic structure of the system is more complicated.

.. image:: Si_w.png
       :height: 400 px

Good convergence is reached for :math:`\omega_\text{lin} \approx \omega_\text{max}/3` and :math:`\Delta\omega` = 0.05 eV.

.. note::

    If more memory is needed, use ``wpar=int`` to parallelize over frequency points. ``int`` needs to be an integer divisor of the available cores.
