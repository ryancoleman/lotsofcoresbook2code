.. _gw exercise:

===========================================
G0W0 calculation of the band gap of silicon
===========================================

In this exercise we will calculate the band gap of silicon.  Go through
the :ref:`bandstructures` tutorial first. In that
tutorial the band structure of silicon is calculated based on the Kohn-Sham
eigenvalues obtained from a DFT calculation. Often one interprets the band
structure as electron addition and removal energies and the difference
between the conduction band minimum and valence band maximum as the band gap.
With the values obtained from DFT this is however not completely correct,
because DFT only finds the ground state electron density when the system is
neutral. When we add or remove an electron we charge the system and the other
electrons redistribute to a new configuration. Thus DFT tends to
underestimate the band gap, since it does not take into account the extra
energy it requires to squeeze an extra electron into the conduction band and
the relieve of energy when we remove an electron from the valence band.

* How many valence electron does silicon have? What are the band indices of
  the valence and conduction bands respectively? Take a look at the band
  structure of silicon - where are the valence band maximum and conduction band
  minimum located and what is the band gap obtained from DFT?

The GW approximation is a method for calculating the charged states of a
system based on a systematic pertubation theory approach. In its simplest
version, and the one so far implemented in GPAW, we can use it to find the
corrections to the DFT band structure than include the missing screening
effect when we add or remove an electron. In this approximation the
quasiparticle energy (electron addition/removal energies) of the state
:math:`(\mathbf{k}, n)` is given by

.. math::
    
    \epsilon^\text{qp}_{n \mathbf{k}} =
    \epsilon_{n \mathbf{k}} + Z_{n \mathbf{k}} \cdot \text{Re}
    \left(\Sigma_{n \mathbf{k}}^{\vphantom{\text{XC}}} +
    \epsilon^{\text{EXX}}_{n \mathbf{k}} -
    V^{\text{XC}}_{n \mathbf{k}} \right).

For more information on the theory of GW, there's a short description in the
section: :ref:`gw_theory`.


------------------------
Ground state calculation
------------------------

As GW is a pertubative approach on top of a ground state, lets start by
calculating the ground state of Silicon. You can create a script yourself or
reuse the one from the band structure calculation. Silicon has in diamond
structure with a lattice constant of 5.431 Å.

.. note:: Currently the implementation of GW in GPAW does not support magnetic systems, so make sure your ground state calculation is performed without spin polarization.

In order to carry out the GW calculation we need the wavefunctions and corresponding energies of both all the occupied and a lot of unoccupied bands. By default, GPAW only calculates the occupied bands, since they are the only ones relevant for the ground state, so you have to specify that you want unoccupied bands as well. This can be done by specifying the parameters ``nbands=XX`` and ``convergence={'bands': 'all'}``, however with the standard iterative eigensolvers it can be hard to converge all the bands. Thus a better way is to first calculate the ground state and the occupied bands and then calculate the unoccupied bands by exact diagonalization of the resulting Hamiltonian. This is done like this:

::

  calc = GPAW(mode=PW(300), ...) # Note we MUST use plane waves

  atoms.set_calculator(calc)
  atoms.get_potential_energy()   # Calculate ground state

  calc.diagonalize_full_hamiltonian() # Does what it says ;)
  calc.write('my_file.gpw', 'all')    # Saves calculator AND wavefunctions

Note that in order to use the ``diagonalize_full_hamiltonian()`` method, the calculator has to use a plane wave basis, which is done by specifying ``mode=PW()`` and correctly importing the PW class: ``from gpaw import PW``. It is a good idea to save the calculator with the ``write()`` method and remember to include the wavefunctions by setting the second argument to ``'all'``.

- Write a script that calculates the ground state of silicon and all of the unoccupied bands and saves the result. Use a plane wave cut-off of say 200 eV and 4x4x4 k-point sampling that is shifted to include the :math:`\Gamma`-point (this can be done by specifying ``kpts={'size': (4, 4, 4), 'gamma': True}``. Why this k-point sampling?

----------------
G0W0 calculation
----------------

To do a GW calculation is very easy. First we must decide which states we actually want to perform the calculation for. For just finding the band gap we can many times just do with the locations of the conduction band minimum and valence band maximum. However the quasiparticle spectrum might be qualitatively different from the DFT spectrum, so its best to do the calculation for all k-points. Here's a script that does this:

.. literalinclude:: Si_g0w0_ppa.py

- Try using the above script to calculate the GW band gap. Is it larger than the one from DFT as we would expect?

- The G0W0 calculator generates a couple of temporary files and output files. Take a look in the output files (ending in ``.txt``) to see what kind of information that is stored.

In the above script we told the GW calculator to use ``ppa=True``, which means that it will use the Plasmon Pole Approximation. In GW we take the full dynamical screening into account and we thus need the dielectric function for the full frequency range. However in the Plasmon Pole Approximation the frequency dependence of the dielectric function is obtained by fitting to a special function that only needs two values to be determined, making the calculation much faster. If we disable the plasmon pole approximation, the dielectric function will be calculated on a frequency grid.

- Try doing the same calculation with ``ppa=False`` and notice the time it takes. Are the results much different from the ones obtained with the Plasmon Pole Approximation?

-----------
Convergence
-----------
As with all other numerical implementations, we have to be careful to check whether the calculations are converged with respect to the various parameters. In GW the most crucial parameters are the density of points that we integrate over: k-points and frequency points, as well as the number of bands included in the band summations and the quality of the wavefunctions, determined by the plane wave cut-off energy.

- By default GPAW chooses the number of bands close to the number of plane waves. It is usually not necesarry to include more bands than the size of the plane wave basis used to expand the wavefunctions. Can you give a simple explanation why?

The fine features of the dielectric function are often averaged out in the integral over the screened potential, which is why the plasmon pole approximation usually performs well.

- Try making the frequency grid denser or coarser by setting the parameter ``domega0`` to something different than its default value, which is ``domega0=0.025``. The calculation time scales linearly with the number of frequency points, so making it half as big doubles the time. Do your results depend a lot on the frequency grid? When can you safely say they are converged?

Next, we need to make sure that we have enough plane waves to properly describe the wavefunctions, by adjusting the plane wave cut-off. This, however, does not come for free; The GW calculation scales quadratically in the number of plane waves and since we set the number of bands to the same number, it actually scales as the **third power** of the number of plane waves!!!

- Try making a couple of calculations where you change the plane wave cut-off energy from say 25 eV to 150 eV. Just use the default frequency grid. If you wish, you can actually read off the number of plane waves used by looking in the generated log file for the screened potential that ends in ``.w.txt``.

Lastly, we also need to make sure that the calculations is converged with respect to the k-point sampling. To do this, one must make new ground state calculations with different k-point samplings to be put into the G0W0 calculator. The calculation of the quasiparticle energy of one state scales quadratically in the number of k-points, but if one want the full band structure there's an extra factor of the number of k-points, so this quickly becomes very heavy!

- Make new groundstate calculations with k-point samplings 4x4x4, 6x6x6 and 8x8x8 and so on and find the DFT band gap. When is it safe to say that the DFT band gap is converged?
- Perform GW calculations (parallelize over minimum four cpus) for the different k-point samplings (4, 6 and 8 k-points only) and compare the gaps. How big is the variation in the gaps compared to the variation in the DFT result? When do you think the GW band gap is converged?
