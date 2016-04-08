BSE for excitonic effects
=========================

Excitonic effects play a fundamental role in determining the optical spectum
of a semiconductor. As soon as an electron from the valence band is promoted
to the conduction band by an external electromagnetic field, it interacts
with the hole left behind in valence band via Coulomb interaction; the
electron-hole pair originated from this process is refferred as exciton.
Eventually the exciton contributes with accessible states inside the band gap
of the semicondutor shifting the onset of the optical absorption and
producing new charecteristic features in the spectrum. Being the exciton a
many body effect involving two particles, it cannot be described by a simple
single-particle picture.  The most refined ab-initio method for calculating
the absorption spectrum is consists in the numerical solution of the Bethe-
Salpeter equation (BSE). In this exercise we show how this can be done using
GPAW. For a brief overview on the BSE and how it is implemented in GPAW take
a look at this page: :ref:`bse`.

In the following we will study the absorption spectrum of Lithium Flouride,
which is well known for having strong excitonic effects. As a starting point
for a BSE calculation we need the single-particle eigenvalues and
eigenfunctions for a large number of bands and k-points. However, since BSE
calculations are extremely heavy in the exercise we use parameters which are
not fully converged in order to keep the calculation doable in a reasonable
amount of time. Nevertheless it is important to remind that before comparing
the results with experiment or other codes, all the convergence parameters
have to be carefully checked. To obtain eigenvalues and eigenfunctions we
first perform a self-consistent calculation for getting the ground-state
density and then we diagonalize the full single-particle Hamiltonian for
getting eigenvalues and eigenenergies also for the excited states. Download
the script :download:`LiF_gs.py`.
 
.. literalinclude:: LiF_gs.py

Read it, try to understand what it does and run it::

    python LiF_gs.py

Before proceding with the BSE we want to calculate the single-particle
spectrum in order to have something to compare with. We do this using the RPA
approximation (as you probably did in the previous exercises).
Download the script :download:`LiF_RPA.py`.

.. literalinclude:: LiF_RPA.py

Once you understand it run it in parallel to save time::

    mpirun -np 8 gpaw-python LiF_RPA.py

The script produces a .csv file containing the imaginary part of the
dielectric function including local field effects. We can take a look to the
absorption spectrum plotting the first and the last columns in the file
``df.csv``.

The absorption spectrum just produced includes only single-particle
transitions and therefore the onset of absorption have to correspond to the
electronic band gap. Try to plot the band structure of LiF and verify that
this is the case.

Now we can go on with the BSE calculation running the script
:download:`LiF_BSE.py`:

.. literalinclude:: LiF_BSE.py

Even if the parameter chosen are not converged, the calculation takes quite
some time; therefore we would better submit the script to nifhleim and run it
on 32 processors typing::

    gpaw-qsub -q small -l nodes=4:ppn=8:xeon8 LiF_BSE.py

The file df.dat finally contains the absorption spectrum including excitonic
effects (if you do not have enough time left you can see the results directly
here: :download:`df.dat`). Try to plot it on top of the single-particle
spectrum, considering that the df.dat file contains the frequencies as first
column and the imaginary part of the dielectric function as third.

Can you tell which are the main differences? Is the onset of the optical
transition still equal to the electronic band gap? What about the single
particle regime, is it anyhow affected by the inclusion of excitonic effects?
How and why?

An important quantity that one could evaluate from the absorption spectrum is
the exciton binding energy (which sets the onset of the optical transitions).
It is defined according to the following expression:

.. math::
    
    E_b = E_{gap}-E_{optical onset}

The experimental value for LiF is ~1.6 eV What can you then say about the
calculation just performed?
