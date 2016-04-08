.. _timepropagation:

======================
Time-propagation TDDFT
======================

Optical photoabsorption spectrum as well as nonlinear effects can be
studied using time propagation TDDFT. This approach
scales better than linear response, but the prefactor is so large that
for small and moderate systems linear response is significantly
faster.


------------
Ground state
------------

To obtain the ground state for TDDFT, one has to just do a standard ground state 
with a larger simulation box. A proper distance from any atom to edge of the 
simulation box is problem dependent, but a minimum reasonable value is around
6 Ångströms and recommended between 8-10 Ång. In TDDFT, one can use larger 
grid spacing than for geometry optimization. For example, if you use h=0.25
for geometry optimization, try h=0.3 for TDDFT. This saves a lot of time. 

A good way to start is to use too small box (vacuum=6.0), too large grid 
spacing (h=0.35), and too large time step (dt=16.0). Then repeat the simulation
with better parameter values and compare. Probably lowest peaks are already 
pretty good, and far beyond the ionization limit, in the continuum, the spectrum 
is not going to converge anyway. The first run takes only fraction of 
the time of the second run.

For a parallel-over-states TDDFT calculation, you must choose the number 
of states so, that these can be distributed equally to processors. For 
example, if you have 79 occupied states and you want to use 8 processes 
in parallelization over states, add one unoccupied state to get 80 states 
in total.


Ground state example::

  # Standard magic
  from ase import Atoms
  from gpaw import GPAW
  
  # Beryllium atom
  atoms = Atoms(symbols='Be', 
                positions=[(0, 0, 0)],
                pbc=False)
  
  # Add 6.0 ang vacuum around the atom
  atoms.center(vacuum=6.0)
  
  # Create GPAW calculator
  calc = GPAW(nbands=1, h=0.3)
  # Attach calculator to atoms
  atoms.set_calculator(calc)
  
  # Calculate the ground state
  energy = atoms.get_potential_energy()
  
  # Save the ground state
  calc.write('be_gs.gpw', 'all')



--------------------------------
Optical photoabsorption spectrum
--------------------------------

Optical photoabsorption spectrum can be obtained by applying a weak 
delta pulse of dipole electric field, and then letting the system evolve
freely while recording the dipole moment. A time-step around 4.0-8.0
attoseconds is reasonable. The total simulation time should be few tens
of femtoseconds depending on the desired resolution.


Example::

  from gpaw.tddft import *
  
  time_step = 8.0                  # 1 attoseconds = 0.041341 autime
  iterations = 2500                # 2500 x 8 as => 20 fs
  kick_strength = [0.0,0.0,1e-3]   # Kick to z-direction
  
  # Read ground state
  td_calc = TDDFT('be_gs.gpw')
  
  # Kick with a delta pulse to z-direction
  td_calc.absorption_kick(kick_strength=kick_strength)
  
  # Propagate, save the time-dependent dipole moment to 'be_dm.dat',
  # and use 'be_td.gpw' as restart file
  td_calc.propagate(time_step, iterations, 'be_dm.dat', 'be_td.gpw')

  # Calculate photoabsorption spectrum and write it to 'be_spectrum_z.dat'
  photoabsorption_spectrum('be_dm.dat', 'be_spectrum_z.dat')

.. note::

  Make sure to number of iterations is divisible by the dump interval
  such that the last iteration will be stored in the restart file.
  Otherwise append td_calc.write('be_td.gpw', mode='all') to the script.

When propagating after an absorption kick has been applied, it is a good
idea to periodically write the time-evolution state to a restart file.
This ensures that you can resume adding data to the dipole moment file
if you experience artificial oscillations in the spectrum because the total
simulation time was too short.

Example::

  from gpaw.tddft import *
  
  time_step = 8.0                  # 1 attoseconds = 0.041341 autime
  iterations = 2500                # 2500 x 8 as => 20 fs

  # Read restart file with result of previous propagation
  td_calc = TDDFT('be_td.gpw')

  # Propagate more, appending the time-dependent dipole moment to the
  # already existing 'be_dm.dat' and use 'be_td2.gpw' as restart file
  td_calc.propagate(time_step, iterations, 'be_dm.dat', 'be_td2.gpw')

  # Recalculate photoabsorption spectrum and write it to 'be_spectrum_z2.dat'
  photoabsorption_spectrum('be_dm.dat', 'be_spectrum_z2.dat')


Typically in experiments, the spherically averaged spectrum is measured.
To obtain this, one must repeat the time-propagation to each Cartesian 
direction and average over the Fourier transformed dipole moments.


--------------------------------
Fourier transformed density
--------------------------------

If one merely wishes to record the time-evolution of the dipole moment and
analyze the resulting spectrum, TDDFT offers little advantage over the much
faster :ref:`LrTDDFT <lrtddft>`. Further, one must bear in mind that only
excitations induced by the absorption kick will show up in the spectrum.

However, propagating a slightly perturbed ground state density may offer
much more structural information, starting with the ability to distinguish
which spectral peaks correspond to which principal directions in a lattice.

Since the dipole moment is generated by displacements in the charge density, 
most strong peaks in the optical photoabsorption spectrum signify nearly
harmonic oscillations herein. Therefore, taking Fourier transforms of the
time-evolution of the density at the resonant frequencies is a great way of
analyzing the spatial extent of the oscillating modes.


The discrete moving-average Fourier transform of the pseudo-electron density
:math:`\tilde{n}(\mathbf{r},t)` is defined:

.. math::

    F_N(\mathbf{r},\omega) = \frac{1}{\sqrt{\pi}} \sum_{j=0}^N \big(
    \tilde{n}(\mathbf{r},t_j)-\overline{n}_N(\mathbf{r})\big)
    \mathrm{e}^{-\textstyle\frac{1}{2}t_j^2\sigma^2}
    \mathrm{e}^{\displaystyle\mathrm{i}\omega t_j} \Delta t_j

, where we allow for variable time-step :math:`\Delta t_j` along the :math:`N`
propagation steps in the time-series :math:`j=0,1,\ldots,N`. With a total
propagation time of :math:`t_N`, the Fourier transforms are taken relative to
the time-average :math:`\overline{n}_N(\mathbf{r})` of the pseudo density:

.. math::

    \overline{n}_N(\mathbf{r}) = \frac{1}{t_{N+1}} \sum_{j=0}^N
    \tilde{n}(\mathbf{r},t_j) \Delta t_j \qquad, t_N = 
    \sum_{j=0}^{N-1}\Delta t_j


Regrettably, having arrived at time :math:`t_N` will not enable us to perform
the above summations because recording :math:`N\sim 10^4` sets of grid data is
completely intractable. Instead, an iterative cumulation scheme is implemented,
which only requires data from one time-step at a time.

XXX more on this later


The class :epydoc:`DensityFourierTransform <gpaw.tddft.fourier>` is used to
calculate and maintain Fourier transforms of the pseudo electron density. It
functions by attaching itself to a TDDFT instance, which in turn notifies
it after each time-step and allows it to update the density Fourier transforms.

.. important::

    An incontestable restriction of the iterative approach is the requirement
    that the frequencies must be given upon initialization (i.e. time zero).
    To avoid wasted effort, getting the peak frequencies right is essential.

It is recommended to use either :ref:`LrTDDFT <lrtddft>` or a somewhat cruder 
time-propagation to estimate which frequencies could be of interest. In the
latter case, applying a weak kick ``[1e-3, 1e-3, 1e-3]`` will probably be
sufficient to excite and detect all the relevant modes in a short time-span.
For quick estimates, using the ``ECN`` propagator and the ``CSCG`` eigensolver
with a tolerance around ``1e-4`` works reasonably well for timesteps of 5-10 as.

.. tip::

    Using a finite width :math:`\sigma` around ``0.1 eV`` will make any
    ballpark figure a much safer bet. Be aware that peaks found using
    :ref:`LrTDDFT <lrtddft>` may shift slightly.


Example::

  from gpaw.tddft import TDDFT
  from gpaw.tddft.fourier import DensityFourierTransform

  time_step = 4.0                  # 1 attoseconds = 0.041341 autime
  iterations = 5000                # 5000 x 4 as => 20 fs
  kick_strength = [0.0,5e-3,0.0]   # Kick to y-direction
  frequencies = [4.26,6.27,13.0, \
                 16.9,18.1,19.9]   # Pre-determined peak frequencies in eV
  sigma = 0.05                     # Width of Gaussian envelope in eV

  # Read ground state
  td_calc = TDDFT('bda_gs.gpw')
  
  # Kick with a delta pulse to y-direction
  td_calc.absorption_kick(kick_strength=kick_strength)

  # Create and attach Fourier transform observer
  obs = DensityFourierTransform(timestep, frequencies, sigma)
  obs.initialize(td_calc)

  # Propagate, save the time-dependent dipole moment to 'bda_dm.dat',
  # (just for comparison) and use 'bda_td.gpw' as restart file
  td_calc.propagate(time_step, iterations, 'bda_dm.dat', 'bda_td.gpw')

  # Save result of the Fourier transformations to a .ftd file
  obs.write('bda_fourier.ftd')

You can now resume adding data to both the dipole moment file and the density
fourier transform if the spectrum is not sufficiently evolved because the total
simulation time was too short.

Example::

  from gpaw.tddft import TDDFT
  from gpaw.tddft.fourier import DensityFourierTransform

  time_step = 4.0                  # 1 attoseconds = 0.041341 autime
  iterations = 5000                # 5000 x 4 as => 20 fs
  frequencies = [4.26,6.27,13.0, \
                 16.9,18.1,19.9]   # Pre-determined peak frequencies in eV
  sigma = 0.05                     # Width of Gaussian envelope in

  # Read restart file with result of previous propagation
  td_calc = TDDFT('bda_td.gpw')

  # Create and attach Fourier transform observer
  obs = DensityFourierTransform(timestep, frequencies, sigma)
  obs.initialize(td_calc)

  # Read previous result of the corresponding Fourier transformations
  obs.read('bda_fourier.ftd')

  # Propagate more, appending the time-dependent dipole moment to the
  # already existing 'bda_dm.dat' and use 'bda_td2.gpw' as restart file
  td_calc.propagate(time_step, iterations, 'bda_dm.dat', 'bda_td2.gpw')

  # Save result of the improved Fourier transformations to an .ftd file
  obs.write('bda_fourier2.ftd') 


--------------------------------
Time propagation
--------------------------------

Since the total CPU time also depends on the number of iterations performed
by the linear solvers in each time-step, smaller time-steps around 2.0-4.0
attoseconds might prove to be faster with the ``ECN`` and ``SICN``
propagators because they have an embedded Euler step in each predictor step:

.. math::

  \tilde{\psi}_n(t+\Delta t) \approx (1 - i \hat{S}^{\;-1}_\mathrm{approx.}(t) \tilde{H}(t) \Delta t)\tilde{\psi}_n(t)

, where :math:`\hat{S}^{\;-1}_\mathrm{approx.}` is an inexpensive operation
which approximates the inverse of the overlap operator :math:`\hat{S}`. See
the :ref:`Developers Guide <overlaps>` for details.


Therefore, as a rule-of-thumb, choose a time-step small enough to minimize the
number of iterations performed by the linear solvers in each time-step, but
large enough to minimize the number of time-steps required to arrive at the
desired total simulation time.


--------------------------------
TDDFT reference manual
--------------------------------

The :class:`~gpaw.tddft.TDDFT` class and keywords:

===================== =============== ============== =====================================
Keyword               Type            Default        Description
===================== =============== ============== =====================================
``ground_state_file`` ``string``                     Name of the ground state file
``td_potential``      ``TDPotential`` ``None``       Time-dependent external potential
``propagator``        ``string``      ``'SICN'``     Time-propagator (``'ECN'``/``'SICN'``/``'SITE'``/``'SIKE'``)
``solver``            ``string``      ``'CSCG'``     Linear equation solver (``'CSCG'``/``'BiCGStab'``)
``tolerance``         ``float``       ``1e-8``       Tolerance for linear solver
===================== =============== ============== =====================================

Keywords for :func:`~gpaw.tddft.TDDFT.absorption_kick`:

================== =============== ================== =====================================
Keyword            Type            Default            Description
================== =============== ================== =====================================
``kick_strength``  ``float[3]``    ``[0,0,1e-3]``     Kick strength
================== =============== ================== =====================================

Keywords for :func:`~gpaw.tddft.TDDFT.propagate`:

====================== =========== =========== ================================================
Keyword                Type        Default     Description
====================== =========== =========== ================================================
``time_step``          ``float``               Time step in attoseconds (``1 autime = 24.188 as``)
``iterations``         ``integer``             Iterations
``dipole_moment_file`` ``string``  ``None``    Name of the dipole moment file
``restart_file``       ``string``  ``None``    Name of the restart file
``dump_interal``       ``integer`` ``500``     How often restart file is written
====================== =========== =========== ================================================

Keywords for :func:`gpaw.tddft.photoabsorption_spectrum`:

====================== ============ ============== ===============================================
Keyword                Type         Default        Description
====================== ============ ============== ===============================================
``dipole_moment_file`` ``string``                  Name of the dipole moment file
``spectrum_file``      ``string``                  Name of the spectrum file
``folding``            ``string``   ``Gauss``      Gaussian folding (or Lorentzian in future)
``width``              ``float``    ``0.2123``     Width of the Gaussian/Lorentzian (in eV)
``e_min``              ``float``    ``0.0``        Lowest energy shown in spectrum (in eV)
``e_max``              ``float``    ``30.0``       Highest energy shown in spectrum (in eV)
``delta_e``            ``float``    ``0.05``       Resolution of energy in spectrum (in eV)
====================== ============ ============== ===============================================

.. autoclass:: gpaw.tddft.TDDFT
   :members:

.. autofunction:: gpaw.tddft.photoabsorption_spectrum
