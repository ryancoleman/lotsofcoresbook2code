.. _bandgab:

==================================================
Calculating band gap using the GLLB-sc functional
==================================================

In this tutorial, we use the GLLB-sc to calculate the band gap of KTaO3 using the 
XC functional GLLB-sc. This functional uses the GLLB response potential to 
replace the PBEsol response potential of the exchange. [GLLB-sc]
This has been shown to improve the band gap description as shown in the figure 
below taken from [Castelli2012].

.. figure:: GLLB-SC_gap.png

A GLLB-sc band gap calculation is performed as given here: 

.. literalinclude:: gllbsc_band_gap.py

Spin-polarized GLLB-SC
======================

Spin-polarized GLLB-SC is currently implemented to svn trunk. However there are some convergence
issues releated to fermi smearing and the reference energy of highest orbital. Also some parts are 
still untested. The code will be improved to robust version soon,
but in the meanwhile please contact mikael.kuisma@tut.fi before using.


-------------

.. [GLLB-sc] M. Kuisma, J. Ojanen, J. Enkovaara, and T. T. Rantala1,
   PHYSICAL REVIEW B 82, 115106 (2010),
   *Kohn-Sham potential with discontinuity for band gap materials*,
   DOI: 10.1103/PhysRevB.82.115106

.. [Castelli2012] Ivano E. Castelli, Thomas Olsen, Soumendu Datta,
   David D. Landis, Søren Dahl, Kristian S. Thygesena
   and Karsten W. Jacobsen,
   Energy Environ. Sci., 2012, 5, 5814,
   *Computational screening of perovskite metal oxides for optimal solar light
   capture*,
   DOI: 10.1039/c1ee02717d
