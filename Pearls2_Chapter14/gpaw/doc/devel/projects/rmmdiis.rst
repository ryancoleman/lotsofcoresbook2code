Improving the RMM-DIIS eigensolver
==================================

:Who:
    Jens Jørgen

Currently, our :ref:`RMM-DIIS eigensolver <RMM-DIIS>` will always take
two steps for each state.  In an attempt to make the eigensolver
faster and more robust, we should investigate the effect of taking a
variable number of steps depending on the change in the eigenstate
error and occupations number as described in [Kresse]_.


.. [Kresse] G. Kresse, J. Furthmüller:
   Phys. Rev. B 54, 11169 - 11186 (1996)
   "Efficient iterative schemes for ab initio total-energy calculations
   using a plane-wave basis set"

Any potential improvements in algorithms need to be evaluated on a large
set of systems. The :ref:`scf_conv_eval` page is dedicated to
presenting those evaluations.
