.. _band_parallelization:

====================
Band parallelization
====================

The orthogonalization can be paralleized over **k**-points, spins,
domains (see :ref:`orthogonalization`), and bands, described below.

Let's say we split the bands in five groups and give each group of
wave functions to one of five processes:

The overlap matrix contains 5x5 blocks.  These are the steps::

  rank:      1           2           3           4           5

         A . . . .   . B . . .   . . C . .   . . . . .   . . . . .
         . . . . .   . A . . .   . . B . .   . . . C .   . . . . .
  S:     . . . . .   . . . . .   . . A . .   . . . B .   . . . . C
         C . . . .   . . . . .   . . . . .   . . . A .   . . . . B
         B . . . .   . C . . .   . . . . .   . . . . .   . . . . A

A. Each process calculates its block in the diagonal and sends a copy
   of its wave functions to the right (rank 5 sends to rank 1).

B. Rank 1 now has the wave functions from rank 5, so it can do the row
   5, column 1 block of `\mathbf{S}`.  Rank 2 can do the row 1, column
   2 block and so on.  Shift wave functions to the right.

C. Rank 1 now has the wave functions from rank 4, so it can do the row
   4, column 1 block of `\mathbf{S}` and so on.

Since `\mathbf{S}` is symmetric, we have all we need::

  A B C . .
  . A B C .
  . . A B C
  C . . A B
  B C . . A

With `B` blocks, we need `(B - 1) / 2` shifts.

Now we can calculate `\mathbf{L}^{-1}` and do the matrix product
`\tilde{\mathbf{\Psi}}_0 \mathbf{L}^{-1}` which requires `B - 1`
shifts of the wave functions.
