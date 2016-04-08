==========
Planewaves
==========

With `N=N_1N_2N_3` grid points: `\br^T=(g_1/N_1,g_2/N_2,g_3/N_3)\mathbf
A`, where `g_c=0,1,...,N_c-1`, we get a plane wave expansion of the wave
funtion as:

.. math::

    \tilde\psi_{k n}(\br) =
    \frac{1}{N} \sum_\bG e^{i(\bG+\bk)\cdot \br}c_{\bk n}(\bG),

where the coefficients are given as:

.. math::

    c_{\bk n}(\bG) = \sum_\br e^{-i(\bG+\bk)\cdot\br}\tilde\psi_{\bk n}(\br)


Exact Exchange
==============

From the pair densities:

.. math::

    \tilde\rho_{\bk_1n_1 \bk_2n_2}(\br) =
    \tilde\psi_{\bk_1n_1}(\br)^* \tilde\psi_{\bk_2n_2}(\br) + ... = \\

    \frac{1}{N^2}
    \sum_{\bG\bG'} e^{i(\bG-\bk_1+\bk_2)\cdot \br}
    c_{\bk_1n_1}(\bG)^* c_{\bk_2n_2}(\bG+\bG') =
    \sum_\bG e^{i(\bG-\bk_1+\bk_2)\cdot \br}C_{\bk_1n_1\bk_2n_2}(\bG),

we get the exact exchange energy:

.. math::

    E_x = -\pi\Omega
    \sum_{\bk_1n_1}
    \sum_{\bk_2n_2}
    f_{\bk_1n_1}f_{\bk_2n_2}
    \sum_\bG
    \frac{|C_{\bk_1n_1\bk_2n_2}(\bG)|^2}{|\bk_1-\bk_2-\bG|^2},

where the weight of a `\bk`-point is included in `f_{\bk n}`.  Let
`E_x'` be defined as the sum above excluding the divergent terms
for `\bk_1=\bk_2` and `\bG=0`.  With

.. math::

    F(\bG)=\frac{e^{-\alpha G^2}}{G^2},

we get (see [#Sorouri]_):

.. math::

    E_x = E_x'
    -\pi\Omega\sum_{\bk_1n_1n_2}f_{\bk_1n_1}f_{\bk_1n_2}
    |C_{\bk_1n_1\bk_1n_2}(0)|^2
    \left(\sum_{\bk_2\bG}F(\bk_1-\bk_2-\bG)-
    \sum_{\bk_2}\sum_{\bG\neq\bk_1-\bk_2}F(\bk_1-\bk_2-\bG)\right).

In the limit of an infinitely dense sampling of the BZ and a not too
small `\alpha`, we get

.. math::

    \sum_{\bk_2\bG}F(\bk_1-\bk_2-\bG)=
    \frac{N_k\Omega}{(2\pi)^3}\int_{\text{BZ}}F(\bk)d\bk=
    \frac{N_k\Omega}{(2\pi)^2}\sqrt{\pi/\alpha},

where `N_k` is the number of `\bk`-points.

Finally:

.. math::

    E_x = E_x'
    -\pi\Omega\sum_{\bk_1n_1n_2}f_{\bk_1n_1}f_{\bk_1n_2}
    |C_{\bk_1n_1\bk_1n_2}(0)|^2\gamma,

where

.. math::

    \gamma = 
    \frac{\Omega}{(2\pi)^2}\sqrt{\pi/\alpha}-
    \sum_{\bk}\sum_{\bG\neq\bk}F(\bk-\bG).

The gradient is:

.. math::

   \frac{\partial E_x}{\partial\tilde\psi_{\bk_1n_1}(\br)}=
   -\pi\Omega\sum_{\bk_2n_2}f_{\bk_1n_1}f_{\bk_2n_2}
   e^{i(\bk_1-\bk_2)\cdot\br}\tilde\psi_{\bk_2n_2}(\br)
   \frac1N\sum_\bG\frac{C_{\bk_1n_1\bk_2n_2}(G)^*}{|\bk_1-\bk_2-\bG|^2}
   e^{-i\bG\cdot\br},

where `1/|\bk_1-\bk_2-\bG|^2` is replaced by `\gamma` for the term where
`\bk_1=\bk_2` and `\bG=0`.
   

.. [#Sorouri] *Accurate and Efficient Method for the Treatment of Exchange in a
   Plane-Wave Basis*,
   A. Sorouri, W.M.C. Foulkes, and N.D.M. Hine,
   J. Chem. Phys. 124, 064105-1 -- 064105-7 (2006)
