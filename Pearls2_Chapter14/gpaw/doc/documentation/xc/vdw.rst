.. _vdw:

========================
vdW-DF and BEEF-vdW
========================


Several vdW-DF [#vdW-DF1a]_ type XC functionals
are implemented selfconsistently
in GPAW, and also the BEEF-vdW [#BEEF-vdW]_ density functional.
The vdW-DF variants include vdW-DF [#vdW-DF1a]_, [#vdW-DF1b]_,
vdW-DF2 [#vdW-DF2]_, optPBE-vdW [#opt-vdW]_, optB88-vdW [#opt-vdW]_,
and C09-vdW [#C09-vdW]_.

The selfconsistent implementation uses the Perez-Soler [#soler]_ FFT
algorithm to evaluate the total energy and potential of the
Rutgers-Chalmers nonlocal correlation, which is originally a
six dimensional integral in real space. However, a non-selfconsistent
method which directly sums up the real-space integral is also available.


Doing a vdW-DF calculation
==================================

The selfconsistent FFT method is highly recommended over the real-space method.
Often, the vdW-DF electron density will be very similar to an ordinary GGA
density, so non-selfconsistent evaluations of a vdW-DF type total energy
using the FFT method is often ok. However, vdW-DF forces obviously require
a selfconsistent potential.

As the examples below illustrate, FFT-based vdW-DF calculations
are most easily done by setting e.g. "xc='vdW-DF'"
in the GPAW calculator object.
However, parameters of the FFT algorithm can be assigned non-default values
by importing the vdW-DF base class.


Selfconsistent vdW-DF calculations
-------------------------------------

>>> from ase import *
>>> from gpaw import GPAW
>>> vdw = 'vdW-DF'
>>> atoms = ...
>>> calc = GPAW(xc=vdw, ...)
>>> atoms.set_calculator(calc)
>>> e = atoms.get_potential_energy()


Perturbative vdW-DF calculations (non-selfconsistent)
--------------------------------------------------------

>>> from gpaw import GPAW
>>> xc = 'vdW-DF'
>>> calc = GPAW('input.gpw')
>>> GGA_energy = calc.get_potential_energy()
>>> vdWDF_diff = calc.get_xc_difference(xc)
>>> vdWDF_energy = GGA_energy + vdWDF_diff

In the above examples, other vdW-DF type functionals can be used
by substituting 'vdW-DF2', 'optPBE-vdW', 'optB88-vdW', or 'C09-vdW'
for 'vdW-DF'.
 

Non-default FFT parameters for vdW-DF calculations
-----------------------------------------------------

A number of parameters determine the spline interpolation of the vdW-DF
nonlocal kernel. These may be assigned non-default values if the vdW-DF base
class is explicitly initialized with new settings.
The example below redefines the number of interpolating cubic splines
(Nalpha) used in a vdW-DF2 calculation.

>>> from ase import *
>>> from gpaw import GPAW
>>> from gpaw.xc.vdw import VDWFunctional
>>> vdw = VDWFunctional('vdW-DF2', Nalpha=24)
>>> atoms = ...
>>> calc = GPAW(xc=vdw, ...)
>>> atoms.set_calculator(calc)
>>> e = atoms.get_potential_energy()


Real-space method vdW-DF
------------------------------------

It is also possible to use the much slower real-space method
for non-selfconsistent evaluations of the nonlocal correlation energy,
which might make sense for (very) small systems.
To use the real-space method one must import a class and set a few parameters:

>>> from gpaw.xc.vdw import VDWFunctional
>>> vdw = VDWFunctional('vdW-DF', fft=False, nspins=1, ncut=0.0005)

where nspins=1 is for spin-paired systems and nspins=2 is used
for spin-polarized calculations. A cutoff, ncut, defines how small a density
must be in order not to be included in the 6D integral.


BEEF-vdW functional
===================

The BEEF-vdW density functional uses the vdW-DF2 nonlocal correlation
energy and potential. It is implemented selfconistently in GPAW.
Furthermore, the BEEF-vdW constructions allows the user to calculate
an estimate of the error to be expected on the quantity calculated
selfconsistently with BEEF-vdW (i.e. an error estimate on relative energies,
not on total energies). This estimate stems from non-selfconsistently
applying an ensemble of XC functionals to BEEF-vdW electron densities.
The ensemble error estimate is then computed from the variance
of the ensemble predictions of the quantity of interest.

Below is an example which calculates the BEEF-vdW binding energy
of molecular H2 (E_bind),
as well as an ensemble estimate of the binding energy error (dE_bind)

>>> from ase import *
>>> from gpaw import GPAW
>>> from ase.dft.bee import BEEFEnsemble
>>> xc = 'BEEF-vdW'
>>> h2 = Atoms('H2',[[0.,0.,0.],[0.,0.,0.75]])
>>> h2.center(vacuum=3)
>>> cell = h2.get_cell()
>>> calc = GPAW(xc=xc)
>>> h2.set_calculator(calc)
>>> e_h2 = h2.get_potential_energy()
>>> ens = BEEF_Ensemble(calc)
>>> de_h2 = ens.get_ensemble_energies()
>>> del h2, calc, ens
>>> h = Atoms('H')
>>> h.set_cell(cell)
>>> h.center()
>>> calc = GPAW(xc=xc)
>>> h.set_calculator(calc)
>>> e_h = h.get_potential_energy()
>>> ens = BEEF_Ensemble(calc)
>>> de_h = ens.get_ensemble_energies()
>>> E_bind = 2*e_h - e_h2
>>> dE_bind = 2*de_h[:] - de_h2[:]
>>> dE_bind = dE_bind.std()


Note that the BEEF_Ensemble module has recently been moved from GPAW
to the ASE package.
The default number of ensemble XC functionals is 2000,
for which well-converged error estimates should be ensured.
Therefore, "de_h2" and "de_h" in the example
are both arrays of 2000 pertubations of a BEEF-vdW total energy.
The syntax "ens.get_ensemble_energies(N)" changes this number to N.
The calculator object input to the BEEF_Ensemble class could of course
stem from a restarted GPAW calculation.

It is very important to calculate
the ensemble statistics correctly. Computing the standard deviation of each
array of total energy pertubations makes little sense, only the standard
deviation of the relative energy pertubations should be used for the
BEEF-vdW ensemble error estimates on a quantity.


.. [#vdW-DF1a] M. Dion, H. Rydberg, E. Schroder, D.C. Langreth, and
   B. I. Lundqvist, Van der Waals density functional for general geometries,
   Physical Review Letters, 92, 246401 (2004)

.. [#BEEF-vdW] J. Wellendorff, K. T. Lundgaard, A. Mogelhoj,
   V. Petzold, D. D. Landis, J. K. Norskov, T. Bligard, and K. W. Jacobsen,
   Physical Review B, 85, 235149 (2012)

.. [#vdW-DF1b] M. Dion, H. Rydberg, E. Schroder, D.C. Langreth, and
   B. I. Lundqvist, Erratum: Van der Waals density functional for
   general geometries, Physical Review Letters, 95, 109902 (2005)

.. [#vdW-DF2] K. Lee, D. E. Murray, L. Kong, B. I. Lundqvist,
   and D. C. Langreth, Higher-accuracy van der Waals density functional,
   Physical Review B, 82, 081101 (2010)

.. [#opt-vdW] J. Klimes, D. R. Bowler, and A. Michaelides,
   Chemical accuracy for the van der Waals density functional,
   Journal of Physics: Condensed Matter, 22, 022201 (2010)

.. [#C09-vdW] V. R. Cooper,
   Van der Waals density functional: An appropriate exchange functional,
   Physical Review B, 81, 161104(R) (2010)
   
.. [#soler] Guillermo Román-Pérez and José M. Soler,
   Efficient Implementation of a van der Waals Density Functional: Application
   to Double-Wall Carbon Nanotubes,
   Physical Review Letters 103, 096102 (2009)
