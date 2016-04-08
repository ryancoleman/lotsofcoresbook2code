.. _orthogonalization:

==================================
Orthogonalizing the wave functions
==================================

Let `\tilde{\Psi}_{nG}` be an element of a wave function matrix
holding the value of `\tilde{\psi}_{n}(\mathbf{r}_G)` (state number
`n` and grid point number `G`).  Then we can write the
:ref:`orthogonality requirement <orthogonality>` like this:

.. math::

   \Delta v
   \tilde{\mathbf{\Psi}}^T \hat{\mathbf{O}} \tilde{\mathbf{\Psi}} =
   \mathbf{1},

where `\Delta v` is the volume per grid point and

.. math::

   \hat{\mathbf{O}} = \mathbf{1} +
   \sum_a \tilde{\mathbf{P}}^a \mathbf{\Delta O}^a
   (\tilde{\mathbf{P}}^a)^T
   \Delta v

is the matrix form of the overlap operator.  This matrix is very
sparse because the projector functions `\tilde{P}^a_{iG} =
\tilde{p}^a_i(\mathbf{r}_G - \mathbf{R}^a)` are localized inside the
augmentation spheres.  The `\Delta O^a_{i_1i_2}` atomic PAW overlap
corrections are small `N_p^a \times N_p^a` matrices (`N_p^a \sim 10`)
defined :ref:`here <overlaps>`.



Gram-Schmidt procedure
======================

The traditional sequential Gram-Schmidt orthogonalization procedure is
not very efficient, so we do some linear algebra to allow us to use
efficient matrix-matrix products.  

Let `\tilde{\mathbf{\Psi}}_0` be the non-orthogonal wave functions.
We calculate the overlap matrix:

.. math::

   \mathbf{S} = 
   \Delta v
   \tilde{\mathbf{\Psi}}_0^T \hat{\mathbf{O}} \tilde{\mathbf{\Psi}}_0,

from the raw overlap `\tilde{\mathbf{\Psi}}_0^T
\tilde{\mathbf{\Psi}}_0` and the projections `(\tilde{\mathbf{P}}^a)^T
\tilde{\mathbf{\Psi}}_0`.

This can be Cholesky factored into `\mathbf{S} = \mathbf{L}^T
\mathbf{L}` and we can get the orthogonalized wave functions as:

.. math::

   \tilde{\mathbf{\Psi}} = \tilde{\mathbf{\Psi}}_0 \mathbf{L}^{-1}.


Parallelization
===============

The orthogonalization can be paralleized over **k**-points, spins,
domains, and bands.


**k**-points and spins
----------------------

Each **k**-point and each spin can be treated separately.


Domains
-------

Each domain will have its contribution to the overlap matrix, and these
will have to be summed up using the domain communicator.  The dense
linear algebra can be performed in a replication fashion on all MPI
tasks using LAPACK or in parallel on a subset of MPI tasks using ScaLAPACK.


Bands
-----

Band parallelization is described at :ref:`Band parallelization <band_parallelization>`.
