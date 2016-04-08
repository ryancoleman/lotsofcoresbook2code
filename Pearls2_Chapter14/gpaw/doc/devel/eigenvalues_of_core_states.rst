.. _eigenvalues_of_core_states:

==========================
Eigenvalues of core states
==========================

Calculating eigenvalues for core states can be useful for XAS, XES and
core-level shift calculations.  The eigenvalue of a core state `k`
with a wave function `\phi_k^a(\mathbf{r})` located on atom number
`a`, can be calculated using this formula:

.. math::

  \epsilon_k = \frac{\partial E}{\partial f_k} =
  \frac{\partial}{\partial f_k}(\tilde{E} - \tilde{E}^a + E^a),

where `f_k` is the occupation of the core state.  When `f_k` is
varied, `Q_L^a` and `n_c^a(r)` will also vary:

.. math::

  \frac{\partial Q_L^a}{\partial f_k} = 
  \int d\mathbf{r} Y_{00}
  [\phi_k^a(\mathbf{r})]^2 \delta_{\ell,0} = Y_{00},

.. math::

  \frac{\partial n_c^a(r)}{\partial f_k} = 
  [\phi_k^a(\mathbf{r})]^2.

Using the PAW expressions for the :ref:`energy
contributions<developersguide_total_energy>`, we get:

.. math::

  \frac{\partial \tilde{E}}{\partial f_k} = 
  Y_{00}
  \int d\mathbf{r}
  \int d\mathbf{r}'
  \frac{\tilde{\rho}(\mathbf{r}')
  \hat{g}_{00}^a(\mathbf{r} - \mathbf{R}^a)}
  {|\mathbf{r} - \mathbf{r}'|}
   =
  Y_{00}
  \int d\mathbf{r}
  \tilde{v}_H(\mathbf{r})
  \hat{g}_{00}^a(\mathbf{r} - \mathbf{R}^a),
 
.. math::

  \frac{\partial \tilde{E}^a}{\partial f_k} = 
  Y_{00}
  \int_{r<r_c^a}d\mathbf{r}
  \int_{r'<r_c^a}d\mathbf{r}'
  \frac{\tilde{\rho}^a(\mathbf{r}')
  \hat{g}_{00}^a(\mathbf{r}) }
  {|\mathbf{r} - \mathbf{r}'|}
 
.. math::

  \frac{\partial E^a}{\partial f_k} = 
  -\frac{1}{2} 
  \int d\mathbf{r}
  \phi_k^a(\mathbf{r})
  \nabla^2 \phi_k^a(\mathbf{r}) +
  \int_{r<r_c^a}d\mathbf{r}
  \int_{r'<r_c^a}d\mathbf{r}'
  \frac{\rho^a(\mathbf{r}')
  [\phi_k^a(\mathbf{r})]^2 }
  {|\mathbf{r} - \mathbf{r}'|} +
  \int_{r<r_c^a}d\mathbf{r}
  \frac{\delta E_{\text{xc}}[n(\mathbf{r})]}
  {\delta n} [\phi_k^a(\mathbf{r})]^2
 
