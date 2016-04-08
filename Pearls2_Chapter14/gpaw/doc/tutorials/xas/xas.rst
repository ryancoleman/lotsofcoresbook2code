==========
XAS theory
==========

Schematic illustration of XAS (from [Nil04]_):

.. figure:: xas_illustration.png
   :width: 250 px

The oscillator strengths are proportional to `|\langle \phi_{1s}|
\mathbf{r} | \psi_n \rangle|^2`, where the one-center expansion of
`\psi_n` for the core-hole atom can be used.


XAS examples
============

First we must create a core hole setup.  This can be done with the
:program:`gpaw-setup` command::

    gpaw-setup -f PBE N --name hch1s --core-hole=1s,0.5

or you can write a small script to do it:

.. literalinclude:: setups.py

Set the location of setups as decribed on :ref:`installationguide_setup_files`.

Spectrum calculation using unoccupied states
============================================

We do a "ground state" calculation with a core hole. Use a lot of
unoccupied states.

.. literalinclude:: run.py

To get the absolute energy scale we do a Delta Kohn-Sham calculation
where we compute the total energy difference between the ground state
and the first core excited state. The excited state should be spin
polarized and to fix the occupation to a spin up core hole and an
electron in the lowest unoccupied spin up state (singlet) we must set
the magnetic moment to one on the atom with the hole and fix the total
magnetic moment with ``occupations=FermiDirac(0.0, fixmagmom=True)``:

.. literalinclude:: dks.py


Plot the spectrum:

.. literalinclude:: plot.py


.. figure:: xas_h2o_spectrum.png
   :width: 400 px


Haydock recursion method
========================

For systems in the condensed phase it is much more efficient to use the Haydock
recursion method to calculate the spectrum, thus avoiding to determine
many unoccupied states. First we do a core hole calculation with
enough k-points to converge the ground state density. Then we compute
the recursion coefficients with a denser k-point mesh to converge the
uncoccupied DOS. A Delta Kohn-Sham calculation can be done for the
gamma point, and the shift is made so that the first unoccupied
eigenvalue at the gamma point ends up at the computed total energy
difference. 


.. literalinclude:: diamond1.py

Compute recursion coefficients:

.. literalinclude:: diamond2.py

Compute the spectrum with the get_spectra method. delta is the HWHM
(should we change it to FWHM???) width of the lorentzian broadening,
and fwhm is the FWHM of the Gaussian broadening. 

::

  sys.setrecursionlimit(10000)
  
  name='diamond333_hch_600_1400a.rec'
  x_start=-20
  x_end=100
  dx=0.01
  x_rec = x_start + npy.arange(0, x_end - x_start ,dx)
  
  r = RecursionMethod(filename=name)
  y = r.get_spectra(x_rec, delta=0.4, fwhm=0.4 )
  y2 = sum(y)
  
  p.plot(x_rec + 278.344673457,y2)
  p.show()


Below the calculated spectrum of Diamond with half and full core holes
are shown along with the experimental spectrum. 

.. figure:: h2o_xas_3.png
   :width: 400 px

.. figure:: h2o_xas_4.png
   :width: 400 px

XES
===

To compute XES, first do a ground state calcualtion with an 0.0 core
hole (an 'xes1s' setup as created above ). The system will not be
charged so the setups can be placed on all atoms one wants to
calculte XES for. Since XES probes the occupied states no unoccupied
states need to be determined. Calculate the spectrum with

::

  xas = XAS(calc, mode='xes', center=n)

Where n is the number of the atom in the atoms object, the center
keyword is only necessary if there are more than one xes setup. The spectrum 
can be shifted by putting the first transition to the fermi level, or to 
compute the total energy diffrence between the core hole state and the state
with a valence hole in HOMO.


Further considerations
======================

For XAS: Gridspacing can be set to the default value. The shape of the
spectrum is quite insensitive to the functional used, the DKS shifts
are more sensitive. The absolute energy position can be shifted so
that the calculated XPS energy matches the expreimental value
[Leetmaa2006]. Use a large box, see convergence with box size for a
water molecule below:

.. literalinclude:: h2o_xas_box1.py

and plot it

.. literalinclude:: h2o_xas_box2.py

.. figure:: h2o_xas_box.png
        :width: 550 px  


.. [Nil04] *Chemical bonding on surfaces probed by X-ray emission
   spectroscopy and density functional theory*, A. Nilsson and
   L. G. M. Pettersson, Surf. Sci. Rep. 55 (2004) 49-167
