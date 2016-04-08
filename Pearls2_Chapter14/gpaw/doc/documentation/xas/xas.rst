.. _xas:

===================================
X-Ray Absorption Spectroscopy (XAS)
===================================

Introduction
============

The pseudo wave-functions are solutions to this generalized eigenvalue problem:

.. math::

  H \tilde{\psi}_n = \epsilon_n S \tilde{\psi}_n.

This can be transformed into a standard eigenvalue problem:

.. math::

  S^{-1/2} H S^{-1/2} \psi_n = \epsilon_n \psi_n,

where `\psi_n = S^{1/2} \tilde{\psi}_n` is an all-electron wave function.


XAS cross section
=================

For the cross section, we need this quantity:

.. math::

  \langle \psi_n | x | \phi^a \rangle =
  \sum_i \langle \tilde{\psi}_n | \tilde{p}_i^a \rangle
  \langle \phi_i^a | x | \phi^a \rangle =
  \langle \tilde{\psi}_n | \tilde{\phi}^a \rangle,

where `\phi^a` is the core state localized on atom `a` and
`\tilde{\phi}^a = \sum_i \langle \phi_i^a | x | \phi^a \rangle
\tilde{p}_i^a`.  Now, the cross section is:

.. math::

  \sum_n |\langle \tilde{\psi}_n | \tilde{\phi}^a \rangle|^2
         \delta(\epsilon_n - E) =
  \sum_n \langle \tilde{\phi}^a | S^{-1/2} | \psi_n \rangle
         \delta(\epsilon_n - E)
         \langle \psi_n | S^{-1/2} | \tilde{\phi}^a \rangle.

By introducing `G(E) = (E - S^{-1/2} H S^{-1/2} + i \gamma)^{-1}`, we
get:

.. math::

  \text{Im}[\langle S^{-1/2} \tilde{\phi}^a | G(E) | S^{-1/2} \tilde{\phi}^a \rangle].
  






Recursion method
================

Instead of working with the `u_i` functions from the Taillefumier
paper, we introduce `w_i=S^{1/2}u_i` which are the actual functions
that we need to find.  We now define `y_i` and `z_i` as:

.. math::

  w_i = S z_i,

.. math::

  y_i = H z_i.

With these definitions, the recursion formula reads:

.. math::

   y_i = a_i w_i + b_{i+1} w_{i+1} + b_i w_{i-1},

where:

.. math::

  a_i = \langle z_i | y_i \rangle,

and

.. math::

  b_i = \langle z_i | y_{i-1} \rangle = \langle z_{i-1} | y_i \rangle.

The `w_i` functions should be normalized as:

.. math::

  \langle w_i | S^{-1} | w_i \rangle = \langle w_i | z_i \rangle = 1,

and the recursion is started with `w_0 \propto \tilde{\phi}^a`.



Inverting the S matrix
======================

The S (or O) operator is defined as:

.. math::

  \hat O = 1 + \sum_a \sum_{i_1 i_2} |\tilde p^a_{i_1}> O^a_{i_1 i_2}< \tilde p^q_{i_2}|
 
Where `O^a_{i_1 i_2} = <\phi ^a_{i_1}| \phi ^a_{i_2}> - <\tilde \phi ^a_{i_1}| \tilde \phi ^a_{i_2}>`

Assume that `\hat O^{-1}` can be written as

.. math::

  \hat O^{-1} = 1 + \sum_a \sum_{i_1 i_2} |\tilde p^a_{i_1}> P^a_{i_1 i_2}< \tilde p^a_{i_2}|

Then according to [P.J. Hasnip et al, Comp. Phys. Comm. 174 (2006) 24-29 ] the coefficients `P^a_{i_1 i_2}` are given by

.. math::

  P^a_{i_1 i_2} = -O^a_{i_1 j} ( 1 + B^a_{kl} O^a_{lm} )^{-1}_{j i_2}       

.. math::

  B^a_{kl} = < \tilde p^a_{k}| \tilde p^a_{l}>

With summation over equal indices (except a). These formulas ignore overlap between projectors on different atoms. The accuracy of the `\hat O^{-1}` operator can be checked for example by doing:

.. math::

  <\tilde \phi_{i_1}| \hat O \hat O^{-1} \hat O |\tilde \phi_{i_2}> - \delta_{i_1 i_2} 

which should be zero for all normalized, orthogonalized `\tilde \phi` 
