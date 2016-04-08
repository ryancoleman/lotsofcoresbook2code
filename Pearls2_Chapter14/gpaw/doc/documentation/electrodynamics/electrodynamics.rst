=========================
Classical electrodynamics
=========================

GPAW can perform classical electrodynamics simulations using quasistatic finite-difference time-domain (QSFDTD) method.
In these calculations you must specify the regions with classically polarizable material, as well as their permittivities
:math:`\epsilon(\mathbf{r}, \omega)`.
The electric field and the polarization charge density are propagated in time under an influence of external perturbation.
The QSFDTD method can be also merged with time-propagation simulation, which yields a hybrid multiscale method.

.. toctree::
   :maxdepth: 2

   qsfdtd
   hybridscheme

