==========
TPSS notes
==========

Kinetic energy density
======================

Inside the augmentation sphere of atom `a` (`r<r_c^a`), we have:

.. math::

  \psi_{\sigma\mathbf{k}n}(\mathbf{r}) =
  \sum_i 
  \phi_i^a(\mathbf{r})
  \langle\tilde{p}_i^a | \tilde{\psi}_{\sigma\mathbf{k}n} \rangle.

The kinetic energy density from the valence electrons will be:

.. math::

  \frac{1}{2}
  \sum_{\mathbf{k}n} f_{\sigma\mathbf{k}n} \sum_{i_1i_2}
  \langle \tilde{\psi}_{\sigma\mathbf{k}n} | \tilde{p}_{i_1}^a \rangle
  \langle \tilde{p}_{i_2}^a | \tilde{\psi}_{\sigma\mathbf{k}n} \rangle
  \mathbf{\nabla}\phi_{i_1}^a \cdot \mathbf{\nabla}\phi_{i_2}^a =
  \frac{1}{2}
  \sum_{i_1i_2} D_{\sigma i_1i_2}^a
  \mathbf{\nabla}\phi_{i_1}^a \cdot \mathbf{\nabla}\phi_{i_2}^a.

Here, we insert `\phi_i^a(\mathbf{r})=Y_L\phi_j^a(r)` and use:

.. math::

  \mathbf{\nabla}\phi_i^a(\mathbf{r}) =
  \mathbf{\nabla}Y_L \phi_j^a(r) +
  Y_L \frac{d \phi_j^a}{dr} \mathbf{r} / r,

to get:

.. math::

  \mathbf{\nabla}\phi_{i_1}^a \cdot \mathbf{\nabla}\phi_{i_2}^a =
  \mathbf{\nabla}Y_{L_1} \cdot \mathbf{\nabla}Y_{L_2} 
  \phi_{j_1}^a(r) \phi_{j_2}^a(r) +
  Y_{L_1} Y_{L_2}
  \frac{d \phi_{j_1}^a}{dr} \frac{d \phi_{j_2}^a}{dr}.

Similar equations hold for the pseudo kinetic energy density.
