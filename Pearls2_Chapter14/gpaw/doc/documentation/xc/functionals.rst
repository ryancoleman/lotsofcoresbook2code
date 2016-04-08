.. _xc_functionals:

====================================
Exchange and correlation functionals
====================================

.. index:: libxc


Libxc
=====

We used the functionals from libxc_.  ...



Calculation of GGA potential
============================


In libxc_ we have (see also "Standard subroutine calls" on ccg_dft_design_)
`\sigma_0=\sigma_{\uparrow\uparrow}`,
`\sigma_1=\sigma_{\uparrow\downarrow}` and
`\sigma_2=\sigma_{\downarrow\downarrow}` with

.. math::

  \sigma_{ij} = \mathbf{\nabla}n_i \cdot \mathbf{\nabla}n_j


.. _libxc: http://www.tddft.org/programs/octopus/wiki/index.php/Libxc

.. _ccg_dft_design: http://www.cse.scitech.ac.uk/ccg/dft/design.html


Uniform 3D grid
===============

We use a finite-difference stencil to calculate the gradients:

.. math::

  \mathbf{\nabla}n_g = \sum_{g'} \mathbf{D}_{gg'} n_{g'}.

The `x`-component of `\mathbf{D}_{gg'}` will be non-zero only when `g`
and `g'` grid points are neighbors in the `x`-direction, where the
values will be `1/(2h)` when `g'` is to the right of `g` and `-1/(2h)`
when `g'` is to the left of `g`.  Similar story for the `y` and `z`
components.

Let's look at the spin-`k` XC potential from the energy expression
`\sum_g\epsilon(\sigma_{ijg})`:

.. math::

  v_{kg} = \sum_{g'} \frac{\partial \epsilon(\sigma_{ijg'})}{\partial n_{kg}}
  = \sum_{g'} 
  \frac{\partial \epsilon(\sigma_{ijg'})}{\partial \sigma_{ijg'}}
  \frac{\partial \sigma_{ijg'}}{\partial n_{kg}}

Using `v_{ijg}=\partial \epsilon(\sigma_{ijg})/\partial \sigma_{ijg}`,
`\mathbf{D}_{gg'}=-\mathbf{D}_{g'g}` and

.. math::

  \frac{\partial \sigma_{ijg'}}{\partial n_{kg}} =
  (\delta_{jk} \mathbf{D}_{g'g} \cdot \mathbf{\nabla}n_{ig'} +
   \delta_{ik} \mathbf{D}_{g'g} \cdot \mathbf{\nabla}n_{jg'}),

we get:

.. math::

  v_{kg} = -\sum_{g'} \mathbf{D}_{gg'} \cdot
  (v_{ijg'} [\delta_{jk} \mathbf{\nabla}n_{ig'} +
             \delta_{ik}  \mathbf{\nabla}n_{jg'}]).


The potentials from the general energy expression
`\sum_g\epsilon(\sigma_{0g}, \sigma_{1g}, \sigma_{2g})` will be:

.. math::

  v_{\uparrow g} = -\sum_{g'} \mathbf{D}_{gg'} \cdot
  (2v_{\uparrow\uparrow g'} \mathbf{\nabla}n_{\uparrow g'} +
   v_{\uparrow\downarrow g'} \mathbf{\nabla}n_{\downarrow g'})

and

.. math::

  v_{\downarrow g} = -\sum_{g'} \mathbf{D}_{gg'} \cdot
  (2v_{\downarrow\downarrow g'} \mathbf{\nabla}n_{\downarrow g'} +
   v_{\uparrow\downarrow g'} \mathbf{\nabla}n_{\uparrow g'}).



PAW correction
==============

Spin-paired case:

.. math::

   \Delta E =
   \sum_g 4 \pi w r_g^2 \Delta r_g
   [\epsilon(n_g, \sigma_g) - \epsilon(\tilde n_g, \tilde\sigma_g)],

where `w` is the weight ...

.. math::

    n_g =
    \sum_{i_ii_2} D_{i_1i_2}
    \phi_{j_1g} Y_{L_1}
    \phi_{j_2g} Y_{L_2}
    + n_c(r_g)
    = \sum_L n_{Lg} Y_L,

where

.. math::

    n_{Lg} =
    \sum_q D_{Lq} n_{qg} + \delta_{L,0} \sqrt{4 \pi} n_c(r_g)

and 

.. math::

   D_{Lq} = \sum_p D_p G_{L_1L_2}^L \delta_{q_p,q} = \sum_p D_p B_{Lpq}.

.. math::

    \mathbf{\nabla} n_g =
    \sum_L Y_L \sum_{g'} D_{gg'} n_{Lg'} \hat{\mathbf{r}} +
    \sum_L \frac{n_{Lg}}{r_g} r \mathbf{\nabla} Y_L =
    a_g \hat{\mathbf{r}} + \mathbf{b}_g / r_g.

Notice that `r \mathbf{\nabla} Y_L` is independent of `r` - just as
`Y_L` is.  From the two contributions, which are orthogonal
(`\hat{\mathbf{r}} \cdot \mathbf{b}_g = 0`), we get

.. math::

    \sigma_g =
    a_g^2 + \mathbf b_g \cdot \mathbf b_g / r_g^2.


.. math::

    \frac{\partial \Delta E}{\partial n_{Lg}} =
    4 \pi w \sum_{g'} r_{g'}^2 \Delta r_{g'}
    \frac{\partial \epsilon}{\partial \sigma_{g'}}
    \frac{\partial \sigma_{g'}}{\partial n_{Lg}}.

Inserting

.. math::

    \frac{\partial \sigma_{g'}}{\partial n_{Lg}} =
    2 a_{g'} Y_L D_{g'g} +
    2 \mathbf b_g \cdot (r \mathbf{\nabla} Y_L) \delta_{gg'} / r_g^2,

we get

.. math::

    \frac{\partial \Delta E}{\partial n_{Lg}} =
    8 \pi w \sum_{g'} r_{g'}^2 \Delta r_{g'}
    \frac{\partial \epsilon}{\partial \sigma_{g'}}
    a_{g'} Y_L D_{g'g} +
    8 \pi w \Delta r_g
    \frac{\partial \epsilon}{\partial \sigma_g}
    \mathbf b_g \cdot (r \mathbf{\nabla} Y_L).


Non-collinear case
------------------

.. math::

    \mathbf{m}_g
    = \sum_L \mathbf{M}_{Lg} Y_L.

.. math::

    n_{\alpha g} = (n_g + \alpha m_g) / 2.

.. math::

    2 \mathbf{\nabla} n_{\alpha g} =
    \mathbf{\nabla} n_g +
    \alpha \sum_L (
    Y_L \sum_{g'} D_{gg'}
    \frac{\mathbf{m}_g \cdot \mathbf{M}_{Lg'}}{m_g} \hat{\mathbf{r}} +
    \frac{\mathbf{m}_g \cdot \mathbf{M}_{Lg}}{m_g r_g}
    r \mathbf{\nabla} Y_L)

.. math::

    =
    (a_g + \alpha c_g) \hat{\mathbf{r}} +
    (\mathbf{b}_g + \alpha \mathbf{d}_g) / r_g.

.. math::

    4 \sigma_{\alpha \beta g} =
    (a_g + \alpha c_g) (a_g + \beta c_g)
    + (\mathbf{b}_g + \alpha \mathbf{d}_g) \cdot
    (\mathbf{b}_g + \beta \mathbf{d}_g) / r_g^2.

.. math::

    \frac{\partial c_g}{\partial \mathbf{M}_{Lg'}} =
    \frac{Y_L}{m_g} (
    D_{gg'} \mathbf{m}_g +
    \delta_{gg'} \mathbf{m}_g' -
    \delta_{gg'} \frac{\mathbf{m}_g \cdot \mathbf{m}_g'}{m_g^2}
    \mathbf{m}_g).

.. math::

    \frac{\partial (\mathbf{d}_g)_\gamma}{\partial \mathbf{M}_{Lg'}} =
    \frac{Y_L \delta_{gg'}}{m_g} (
    \mathbf{m}_g r \nabla_\gamma Y_L +
    \sum_{L'} \mathbf{M}_{L'g} r \nabla_\gamma Y_{L'} -
    \frac{\mathbf{m}_g}{m_g^2}
    \sum_{L'} \mathbf{m}_g \cdot \mathbf{M}_{L'g} r \nabla_\gamma
    Y_{L'}).
