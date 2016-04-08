Atomic PAW setups
=================

.. _setup_matrix_elements_nabla:

Calculating matrix elements of nabla
------------------------------------

This integral is needed for LrTDDFT and response function related
quantities:

.. math::

  \langle\phi_i|\mathbf\nabla|\phi_{i'}\rangle -
  \langle\tilde\phi_i|\mathbf\nabla|\tilde\phi_{i'}\rangle,

where `|\phi_i\rangle = \phi_i(\mathbf r) = \phi_j(r)Y_{\ell
m}(\hat{\mathbf r})`, and `|\tilde\phi_i\rangle = \tilde\phi_i(\mathbf
r) = \tilde\phi_j(r)Y_{\ell m}(\hat{\mathbf r})`.

.. math::

  \langle\phi_i|\mathbf\nabla|\phi_{i'}\rangle =
  \langle\phi_i|\frac{\partial}{\partial r}(\phi_{j'}/r^{\ell'})
  \frac{\partial r}{\partial \mathbf r}
  r^{\ell'}Y_{\ell'm'}\rangle +
  \langle\phi_i|\frac{\phi_{j'}}{r^{\ell'}}
  \mathbf\nabla(r^{\ell'}Y_{\ell'm'})\rangle.

Since we use real-valued spherical harmonics, we have:

.. math::

  \frac{\partial r}{\partial \mathbf r}=
  \hat{\mathbf r}=(x/r,y/r,z/r)=
  \sqrt{\frac{4\pi}{3}}(Y_{1m_x},Y_{1m_y},Y_{1m_z}).

Splitting the integral in radial and angular parts, we get:

.. math::

  \langle\phi_i|\frac{\partial}{\partial x}|\phi_{i'}\rangle =
  \sqrt{\frac{4\pi}{3}}
  \int r^2dr
  \phi_j\frac{\partial}{\partial r}(\phi_{j'}/r^{\ell'})r^{\ell'}
  G_{1m_x,\ell'm'}^{\ell m} +
  \int r^2dr
  \phi_j\phi_{j'}/r
  \int d\hat{\mathbf r}
  Y_{\ell m}r^{1-\ell'}\frac{\partial}{\partial x}
  (r^{\ell'}Y_{\ell'm'}),

where `G_{\ell m,\ell'm'}^{\ell''m''}` are Gaunt coefficents and the
last angular integral has been tabulated as ``Y_LLv`` in the
:svn:`~gpaw/gaunt.py` module.


More stuff
----------

.. autoclass:: gpaw.setup.Setup
   :members:

.. autoclass:: gpaw.setup.Setups
   :members:
