.. _densitymix:

==============
Density Mixing
==============


Pulay Mixing
------------

The density is updated using Pulay-mixing [#Pulay1980]_, [#Kresse1996]_.

Pulay mixing (or direct inversion of the iterative subspace (DIIS))
attempts to find a good approximation of the final solution as a
linear combination of a set of trial vectors `\{n^i\}` generated during
an iterative solution of a problem. If the error associated with a
given solution is given as `\{R^i\}` then Pulay mixing assumes that
the error of a linear combination of the trail vectors is given as the
same linear combination of errors

.. math::

  n_{i+1}=\sum \alpha_i n_i \quad,\quad R_{i+1}=\sum \alpha_i R_i

The norm `R^{i+1}` is thus given as 

.. math::

  \langle R_{i+1}|R_{i+1}\rangle=\bar{\alpha}^T \bar{\bar{R}}\bar{\alpha}

where elements of the matrix is given as `\bar{\bar{R}}_{ij}=\langle
R_{i}|R_{j}\rangle`. The norm can thus be minimized by solving

.. math::

  \frac{\delta \langle R_{i+1}|R_{i+1}\rangle}{\delta
  \bar{\alpha}^T}=2 \bar{\bar{R}}\bar{\alpha}=0

In density mixing the error of a given input density is given as

.. math::

  R_i = n_i^{out}[n_i^{in}]-n_i^{in}

The original Pulay mixing only uses `n_i^{out}` to calculate the
errors and thereby the mixing parameters. To more efficiently cover
solution space it can be an advantage to include them with a certain
weight, given as the input parameter `\beta`.

.. math::

  n_{i+1}^{in}=\sum \alpha_i (n_i^{in}+\beta R_i)


Special Metric
--------------

Convergence is improved by an optimized metric `\hat{M}` for
calculation of scalar products in the mixing scheme, `\langle A | B
\rangle _s = \langle A | \hat{M} | B \rangle`, where `\langle \rangle
_s` is the scalar product with the special metric and `\langle
\rangle` is the usual scalar product.  The metric is based on the
rationale that contributions for small wave vectors are more important
than contributions for large wave vectors [#Kresse1996]_.  Using a
metric that weighs short wave density changes more than long wave
changes can reduce charge sloshing significantly.

It has been found [#Kresse1996]_ that the metric

.. math::

  \hat{M} = \sum_q | q \rangle f_q \langle q |, \quad f_q =
  1 + \frac{w}{q^2}

is particularly useful (`w` is a suitably chosen weight).

This is easy to apply in plane wave codes, as it is local in reciprocal space.
Expressed in real space, this metric is

.. math::

  \hat{M} = \sum_{R R'} | R \rangle f(R' - R) \langle R' |, \quad f(R) =
  \sum_q f_q e^{i q R}

As this is fully nonlocal in real space, it would be very costly to apply.
Instead we use a semilocal stencil with only three nearest neighbors:

.. math::

  f(R) = \begin{cases}
  1 + w/8 & R = 0 \\
  w / 16 & R = \text{nearest neighbor dist.} \\
  w / 32 & R = \text{2nd nearest neighbor dist.} \\
  w / 64 & R = \text{3rd nearest neighbor dist.} \\
  0 & \text{otherwise}
  \end{cases}

which corresponds to the reciprocal space metric

.. math::

  f_q = 1 + \frac{w}{8} (1 + \cos q_x + \cos q_y + \cos q_z +
  \cos q_x \cos q_y + \cos q_y \cos q_z + \cos q_x \cos q_z +
  \cos q_x \cos q_y \cos q_z)

With the nice property that it is a monotonously decaying function
from `f_q = w + 1` at `q = 0` to `f_q = 1` anywhere at the zone
boundary in reciprocal space.

A comparison of the two metrics is displayed in the figure below

.. image:: metric.png
  :align: center


Specifying a Mixing Scheme in GPAW
----------------------------------

Specifying the mixing scheme and metric is done using the ``mixer``
keyword of the GPAW calculator::

  from gpaw import GPAW, Mixer
  calc = GPAW(mixer=Mixer(beta=0.05, nmaxold=5, weight=50.0))

which is the recommended value if the default fails to converge.

The class ``Mixer`` indicates one of the possible mixing schemes.  The
Pulay mixing can be based on:

1. The spin densities seperately, ``Mixer`` (This will *not* work for
   a spinpolarized system, unless the magnetic moment is fixed)
2. The total density, ``MixerSum2``
3. Spin channels seperately for the density matrices, and the summed
   channels for the pseudo electron density, ``MixerSum``
4. The total density and magnetization densities seperately, ``MixerDif``

Where the magnetization density is the difference between the two spin
densities.

All mixer classes takes the arguments ``(beta=0.25, nmaxold=3,
weight=50.0)``. In addition, the ``MixerDif`` also takes
the arguments ``(beta_m=0.7, nmaxold_m=2,
weight_m=10.0)`` which is the corresponding mixing parameters for the
magnetization density.

Here ``beta`` is the linear mixing coefficient, ``nmaxold`` is the
number of old densities used, and ``weight`` is the
weight used by the metric, if any.

MixerDif seems to be a good choice for spin polarized
molecules. MixerSum is sometimes better for bulk systems.

The Mixer and MixerSum classes
------------------------------

.. autoclass:: gpaw.mixer.Mixer
   :members:
   :inherited-members:

.. autoclass:: gpaw.mixer.MixerSum
   :members:
   :inherited-members:

References
----------

.. [#Pulay1980] Pulay, Chem. Phys. Let. **73**, 393 (1980)
.. [#Kresse1996] Kresse, Phys. Rev. B **54**, 11169 (1996)
