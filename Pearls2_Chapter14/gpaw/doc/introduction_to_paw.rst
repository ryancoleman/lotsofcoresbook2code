.. _introduction_to_paw:

===================
Introduction to PAW
===================

A simple example
================

We look at the `2\sigma`\ * orbital of a CO molecule: |ts|

.. |ts| image:: _static/2sigma.png

The main quantity in the PAW method is the pseudo wave-function (blue
crosses) defined in all of the simulation box:

.. math::

  \tilde{\psi}(\mathbf{r}) =  \tilde{\psi}(ih, jh, kh),

where `h` is the grid spacing and `(i, j, k)` are the indices of the
grid points.

.. figure:: _static/co_wavefunctions.png

In order to get the all-electron wave function, we add and subtract
one-center expansions of the all-electron (thick lines) and pseudo
wave-functions (thin lines):

.. math::

  \tilde{\psi}^a(\mathbf{r}) =  \sum_i C_i^a \tilde{\phi}_i^a(\mathbf{r})

.. math::

  \psi^a(\mathbf{r}) =  \sum_i C_i^a \phi_i^a(\mathbf{r}),

where `a` is C or O and `\phi_i` and `\tilde{\phi}_i` are atom
centered basis functions formed as radial functions on logarithmic
radial grid multiplied by spherical harmonics.

The expansion coefficients are given as:

.. math::

  C_i^a = \int d\mathbf{r} \tilde{p}^a_i(\mathbf{r} - \mathbf{R}^a)
  \tilde{\psi}(\mathbf{r}).


Approximations
==============

* Frozen core orbitals.
* Truncated angular momentum expansion of compensation charges.
* Finite number of basis functions and projector functions.


More information on PAW
=======================

You can find additional information on the :ref:`literature` page, or
by reading the `paw note`_

.. _paw note: paw_note.pdf



Script
======

.. literalinclude:: co_wavefunctions.py

