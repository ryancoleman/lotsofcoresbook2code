Symmetry
========

Let `\mathbf A^T=(\mathbf a_0,\mathbf a_1, \mathbf a_2)`, where
`\mathbf a_0`, `\mathbf a_1` and `\mathbf a_2` are the lattice vectors
of the unit cell.

.. note::

  `(\mathbf a_c)_v=\mathbf A_{cv}` is stored in ``gd.cell_cv[c, v]``
  in units of Bohr and in ``atoms.cell[c, v]`` in Å units.

The relation between scaled positions `\mathbf s` and xyz-positions
`\mathbf r` is `\mathbf r=\mathbf A^T\mathbf s`.

A crystal has a set of symmetry operations (``symmetry.op_scc``) in
the form of matrices `\mathbf U` so that the lattice vectors are
transformed to `\mathbf A'=\mathbf U\mathbf A` and `\mathbf r` is
transformed to `\mathbf r'` as:

.. math::

   \mathbf r'=
   \mathbf A'^T\mathbf s=
   \mathbf A^T\mathbf U^T\mathbf s=
   \mathbf A^T\mathbf U^T\mathbf A^{-T}\mathbf r=
   \mathbf M\mathbf r,

where `\mathbf M=\mathbf A^T\mathbf U^T\mathbf A^{-T}`.

If we want to express `\mathbf r'` in terms of the original lattice
vectors (`\mathbf r'=\mathbf A^T\mathbf s'`), we get:

.. math::

   \mathbf s' = \mathbf U^T\mathbf s.

.. note::

   The `\mathbf U` matrices contain only the integers -1, 0 and 1.
   Also note, that if `\mathbf U` is a symmetry operation, then
   `\mathbf U^{-1}` is too.

Let `\tilde\psi_{\mathbf k}(\mathbf r)` be a Bloch wave function and
`\mathbf R` any Bravais lattice vector:

.. math::

   \tilde\psi_{\mathbf k}(\mathbf r+\mathbf R)=
   e^{i\mathbf k^T\mathbf R}\tilde\psi_{\mathbf k}(\mathbf r).

Transforming `\tilde\psi_{\mathbf k}` with our symmetry operation, we
get `\tilde\psi'_{\mathbf k'}(\mathbf r)=\tilde\psi_{\mathbf
k}(\mathbf M\mathbf r)` and:

.. math::

   \tilde\psi'_{\mathbf k'}(\mathbf r+\mathbf R)=
   \tilde\psi_{\mathbf k}(\mathbf M\mathbf r+\mathbf M\mathbf R)=
   e^{i\mathbf k^T\mathbf M\mathbf R}
   \tilde\psi_{\mathbf k}(\mathbf M\mathbf r)=
   e^{i\mathbf k^T\mathbf M\mathbf R}
   \tilde\psi'_{\mathbf k'}(\mathbf r).

From this equation it is seen that `\mathbf k'=\mathbf M^T\mathbf k`.
In terms of scaled k-points `\mathbf q`, where:

.. math:: \mathbf k=2\pi\mathbf A^{-1}\mathbf q,

we get `\mathbf q'=\mathbf U\mathbf q`.


Besides cystal symmetry, there is also time reversal symmetry for all 
systems with no magnetic field. The wavefunction for `{\mathbf k}` 
and `{-\mathbf k}` is related as:

.. math::
   
   \tilde\psi_{-\mathbf k}(\mathbf r) = \tilde\psi^{\ast}_{\mathbf k}(\mathbf r)

If in addition the crystal has inversion symmetry, then the wavefunction should 
satisfy: 

.. math::

   \tilde\psi_{\mathbf k}(\mathbf r) = \tilde\psi_{-\mathbf k}(-\mathbf r) 
   =  \tilde\psi^{\ast}_{\mathbf k}(-\mathbf r)

.. note::

   Time reversal symmetry operation is not included in ``symmetry.op_scc``.

Details of the symmetry object
------------------------------

.. autoclass:: gpaw.symmetry.Symmetry
   :members:
