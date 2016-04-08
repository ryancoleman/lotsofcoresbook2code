.. _RMM-DIIS:

================================
The residual minimization method
================================

Algorithm
---------

1) Initial guess for wave functions (`\tilde{\psi}_n`).

2) Orthogonalize wavefunctions (make sure `\langle
   \tilde{\psi}_n | \hat{S} | \tilde{\psi}_m \rangle = \delta_{nm}`).

3) Calculate density (`\tilde{n}`, `D_{ij}^a`).

4) Calculate potential (`\tilde{v}`, `\Delta H_{ij}^a`).

5) Apply hamiltonian (`\hat{H}\tilde{\psi}_n`).

6) Subspace diagonalization (rotate `\tilde{\psi}_n` so that `\langle
   \tilde{\psi}_n | \hat{H} | \tilde{\psi}_m \rangle = \delta_{nm} \epsilon_n`).

7) Calculate residuals (`R_n = \hat{H}\tilde{\psi}_n - \epsilon_n
   \hat{S}\tilde{\psi}_n`).

8) Improve wave functions using the RMM-DIIS algorithm (see below).

9) Back to (2).



RMM-DIIS step
-------------

For each wave function we calculate the residual:

.. math::

 R_n = (\hat{H} - \epsilon_n \hat{S}) \tilde{\psi}_n

New improved wave function: `\tilde{\psi}_n' = \tilde{\psi}_n +
\lambda \hat{P} R_n`, where `\hat{P}` is a preconditioner_.  Find step
length `\lambda` by minimizing the norm of:

.. math::

 R_n' = (\hat{H} - \epsilon_n \hat{S}) \tilde{\psi}_n'

Since we already have `R_n'`, we might as well use it to take an extra
step (with the same step length as for the first step):

.. math::

  \tilde{\psi}_n \leftarrow \tilde{\psi}_n' + \lambda \hat{P} R_n'
  = \tilde{\psi}_n +
  \lambda \hat{P} R_n + \lambda \hat{P} R_n'


See [Kresse96]_ for details.



.. _preconditioner:

Preconditioning
---------------

.. hhhh

   image:: images/preconditioning.png
   :width: 3cm
   :align: center

The ideal preconditioner would be:

.. math::

 \hat{P} = -(\hat{H} - \epsilon_n \hat{S})^{-1}.

For the short wavelength parts of the residuals, `\hat{H} - \epsilon_n
\hat{S}` will be dominated by the kinetic energy operator, so we have
approximately `\hat{P} \simeq -\hat{T}^{-1}`.

We calculate preconditioned residuals (`\tilde{R}_n = \hat{P} R_n`) by
solving `\hat{T} \tilde{R}_n = -R_n` or equivalently

.. math::

  \frac{1}{2} \nabla^2 \tilde{R}_n = R_n

approximately using multigrid techniques as described in [Briggs95]_.




References
----------

.. [Kresse96] G. Kresse, J. Furthmüller:
   Phys. Rev. B 54, 11169 - 11186 (1996)
   "Efficient iterative schemes for ab initio total-energy calculations
   using a plane-wave basis set"

.. [Briggs95] E. L. Briggs, D. J. Sullivan and J. Bernholc:
   Phys. Rev. B 52, R5471 (1995),
   "Large Scale Electronic Structure Calculations with Multigrid
   Acceleration"
