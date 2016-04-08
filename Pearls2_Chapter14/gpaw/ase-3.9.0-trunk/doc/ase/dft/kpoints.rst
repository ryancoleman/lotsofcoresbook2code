.. module:: ase.dft.kpoints
   :synopsis: Brillouin zone sampling

=======================
Brillouin zone sampling
=======================

The **k**-points are always given relative to the basis vectors of the
reciprocal unit cell.


Monkhorst-Pack
--------------

.. autofunction:: monkhorst_pack

The k-points are given as [MonkhorstPack]_:

.. math::

    \sum_{i=1,2,3} \frac{2n_i -N_i - 1}{2N_i} \mathbf{b}_i,

where `n_i=1,2,...,N_i`, ``size`` = `(N_1, N_2, N_3)` and the
`\mathbf{b}_i`'s are reciprocal lattice vectors.

.. autofunction:: get_monkhorst_pack_size_and_offset

Example:

>>> from ase.dft.kpoints import *
>>> monkhorst_pack((4, 1, 1))
array([[-0.375,  0.   ,  0.   ],
       [-0.125,  0.   ,  0.   ],
       [ 0.125,  0.   ,  0.   ],
       [ 0.375,  0.   ,  0.   ]])
>>> get_monkhorst_pack_size_and_offset([[0, 0, 0]])
(array([1, 1, 1]), array([ 0.,  0.,  0.]))


.. [MonkhorstPack]
    Hendrik J. Monkhorst and James D. Pack:
    *Special points for Brillouin-zone integrations*,
    Phys. Rev. B 13, 5188–5192 (1976) 


Chadi-Cohen
-----------

Predefined sets of **k**-points:

.. data:: cc6_1x1
.. data:: cc12_2x3
.. data:: cc18_sq3xsq3
.. data:: cc18_1x1
.. data:: cc54_sq3xsq3
.. data:: cc54_1x1
.. data:: cc162_sq3xsq3
.. data:: cc162_1x1


Naming convention: ``cc18_sq3xsq3`` is 18 **k**-points for a
sq(3)xsq(3) cell.

Try this:

>>> import numpy as np
>>> import pylab as plt
>>> from ase.dft.kpoints import cc162_1x1
>>> B = [(1, 0, 0), (-0.5, 3**0.5 / 2, 0), (0, 0, 1)]
>>> k = np.dot(cc162_1x1, B)
>>> plt.plot(k[:, 0], k[:, 1], 'o')
[<matplotlib.lines.Line2D object at 0x9b61dcc>]
>>> p.show()

.. image:: cc.png


Special points in the Brillouin zone
------------------------------------

You can find the special points in the Brillouin zone:

>>> from ase.lattice import bulk
>>> from ase.dft.kpoints import ibz_points, get_bandpath
>>> si = bulk('Si', 'diamond', a=5.459)
>>> points = ibz_points['fcc']
>>> G = points['Gamma']
>>> X = points['X']
>>> W = points['W']
>>> K = points['K']
>>> L = points['L']
>>> kpts, x, X = get_bandpath([W, L, G, X, W, K], si.cell)
>>> print len(kpts), len(x), len(X)
50 50 6

.. autofunction:: get_bandpath
