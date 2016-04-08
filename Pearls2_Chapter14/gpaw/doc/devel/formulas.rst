========
Formulas
========


Coulomb
=======

.. math::

    \frac{1}{|\br-\br'|} =
    \sum_\ell \sum_{m=-\ell}^\ell
    \frac{4\pi}{2\ell+1}
    \frac{r_<^\ell}{r_>^{\ell+1}}
    Y_{\ell m}^*(\hat\br) Y_{\ell m}(\hat\br')

or

.. math::

    \frac{1}{r} = \int \frac{d\mathbf{G}}{(2\pi)^3}\frac{4\pi}{G^2}
    e^{i\mathbf{G}\cdot\br}.


Fourier transforms
==================

The Fourier transform of a radial function multiplied by a spherical
harmonic is:

.. math::

    f(G)Y_{\ell m}(\hat G) =
    \int d\br e^{i\mathbf{G}\cdot\br} f(r)Y_{\ell m}(\br),

where
    
.. math::

    f(G) = 4\pi i^\ell \int_0^\infty r^2 dr j_\ell(Gr) f(r).

.. note::

    .. math::

        e^{i \mathbf{G} \cdot \br} =
        4 \pi \sum_{\ell m} i^\ell j_\ell(Gr) Y_{\ell m}(\hat{\br})
        Y_{lm}(\hat{\mathbf{G}}).

The `spherical Bessel function`_ is defined as:

.. math::

    j_\ell(x) =
    \text{Re}\{
    \frac{e^{ix}}{x} \sum_{n=0}^\ell
    \frac{(-i)^{\ell+1-n}}{n!(2x)^n}
    \frac{(\ell+n)!}{(\ell-n)!}
    \}.

This is implemented in this function:

.. autofunction:: gpaw.atom.radialgd.fsbt

.. _spherical Bessel function:
    http://en.wikipedia.org/wiki/Bessel_function
    #Spherical_Bessel_functions:_jn.2C_yn


Gaussians
=========

.. math:: n(r) = (\alpha/\pi)^{3/2} e^{-\alpha r^2},

.. math:: \int_0^\infty 4\pi r^2 dr n(r) = 1

Its Fourier transform is:

.. math::

    n(k) = \int d\br e^{i\mathbf{k}\cdot\br} n(r) =
    \int_0^\infty 4\pi r^2 dr \frac{\sin(kr)}{kr} n(r) =
    e^{-k^2/(4a)}.

With `\nabla^2 v=-4\pi n`, we get the potential:

.. math:: v(r) = \frac{\text{erf}(\sqrt\alpha r)}{r},

and the energy:

.. math::

    \frac12 \int_0^\infty 4\pi r^2 dr n(r) v(r) =
    \sqrt{\frac{\alpha}{2\pi}}.

Note: `\text{erf}(x) \simeq x\sqrt{4/\pi}` for small `x`.


Shape functions
---------------

GPAW uses Gaussians as shape functions for the PAW compensation charges:

.. math::

    g_{\ell m}(\br) =
    \frac{\alpha^{\ell + 3 / 2} \ell ! 2^{2\ell + 2}}
    {\sqrt{\pi} (2\ell + 1) !}
    e^{-\alpha r^2}
    Y_{\ell m}(\hat{\br}).

They are normalized as:

.. math::

    \int d \br g_{\ell m}(\br) Y_{\ell m}(\hat{\br}) r^\ell = 1.


Hydrogen
========

The 1s orbital:

.. math:: \psi_{\text{1s}}(r) = 2Y_{00} e^{-r},

and the density is:

.. math:: n(r) = |\psi_{\text{1s}}(r)|^2 = e^{-2r}/\pi.


Radial Schrödinger equation
===========================

With `\psi_{n\ell m}(\br) = u(r) / r Y_{\ell m}(\hat\br)`, we have the
radial Schrödinger equation:

.. math::

   -\frac12 \frac{d^2u}{dr^2} + \frac{\ell(\ell + 1)}{2r^2} u + v u
   = \epsilon u.

We want to solve this equation on a non-equidistant radial grid with
`r_g=r(g)` for `g=0,1,...`.  Inserting `u(r) = a(g) r^{\ell+1}`, we
get:

.. math::

   \frac{d^2 a}{dg^2} (\frac{dg}{dr})^2 r^2 +
   \frac{da}{dg}(r^2 \frac{d^2g}{dr^2} + 2 (\ell+1) r \frac{dg}{dr}) -
   2 r^2 (v - \epsilon) a = 0.


Including Scalar-relativistic corrections
-----------------------------------------

The scalar-relativistic equation is:

.. math::

   -\frac{1}{2 M} \frac{d^2u}{dr^2} + \frac{\ell(\ell + 1)}{2Mr^2} u -
   \frac{1}{(2Mc)^2}\frac{dv}{dr}(\frac{du}{dr}-\frac{u}{r}) + v u
   = \epsilon u.

where the relativistic mass is:

.. math::

   M = 1 - \frac{1}{2c^2} (v - \epsilon).

With `u(r) = a(g) r^\alpha`, `\kappa = (dv/dr)/(2Mc^2)` and

.. math::
    
    \alpha = \sqrt{\ell^2 + \ell + 1 -(Z/c)^2},
    
we get:

.. math::

   \frac{d^2 a}{dg^2} (\frac{dg}{dr})^2 r^2 +
   \frac{da}{dg}(r^2 \kappa \frac{dg}{dr} + r^2 \frac{d^2g}{dr^2} +
   2 \alpha r \frac{dg}{dr}) +
   [2 M r^2 (\epsilon - v) +
   \alpha (\alpha - 1) - \ell (\ell + 1)
    + \kappa (\alpha - 1) r] a = 0.
