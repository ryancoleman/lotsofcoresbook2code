from __future__ import print_function, division
from math import sin, cos, pi
import numpy as np
from gpaw.fd_operators import GUCLaplace as Laplace
from gpaw.fd_operators import Gradient
from gpaw.grid_descriptor import GridDescriptor
from gpaw.mpi import size

cells = [
    ('distorted hexagonal', 4,
     [(1, 0, 0),
      (1.02 * cos(pi / 3 - 0.02), 1.02 * sin(pi / 3 - 0.02), 0),
      (0, 0, 1.0)]),
    ('hexagonal', 4,
     [(1, 0, 0),
      (0.5, 3**0.5 / 2, 0),
      (0, 0, 1.1)]),
    ('fcc', 6,
     [(0, 1, 1),
      (1, 0, 1),
      (1, 1, 0)]),
    ('fcc-alternative', 6,
     [(1, 0, 0),
      (0.5, 3**0.5 / 2, 0),
      (0.5, 3**0.5 / 6, (2 / 3)**0.5)]),
    ('bcc', 4,
     [(-1, 1, 1),
      (1, -1, 1),
      (1, 1, -1)]),
    ('sc', 3,
     [1.1, 1.02, 1.03]),
    ('distorted sc', 6,
     [(1, 0, 0),
      (0.01, 1, 0),
      (0, 0.02, 1)]),
    ('rocksalt', 6,
     [(2 * np.sqrt(1 / 3), np.sqrt(1 / 8), -np.sqrt(1 / 24)),
      (2 * np.sqrt(1 / 3), -np.sqrt(1 / 8), -np.sqrt(1 / 24)),
      (2 * np.sqrt(1 / 3), 0, np.sqrt(1 / 6))]),
    ('nasty', 6,
     [(1, 0, 0),
      (0.0001, 1.03, 0),
      (0.0001, 0.0001, 1.0)]),
    ('Mike', 6,
     5 * np.array([(5.565 / 28, 0, 0),
                   (0.0001 / 28, 5.565 / 28, 0),
                   (0.0001 / 24, 0.0001 / 24, 4.684 / 24)])),
    ('MnO', 6,
     [(1, 0.5, 0.5), (0.5, 1, 0.5), (0.5, 0.5, 1)]),
    ('Jelver', 6,
     [[6.658433626216136, 0.1711724130951983, -0.04300038284455887],
      [3.4774564712755938, 5.843379292501022, 0.01599293966594096],
      [-0.10777038306906983 * 0.43, 0.10850460815311265 * 0.43,
       15.26098014321118 * 0.43]])]


if size == 1:
    for name, D, cell in cells:
        if name == 'Jelver':
            # Strange one!
            continue
            
        print('------------------')
        print(name, D)
        print(cell[0])
        print(cell[1])
        print(cell[2])
        for n in range(1, 6):
            N = 2 * n + 2
            gd = GridDescriptor((N, N, N), cell)
            b_g = gd.zeros()
            r_gv = gd.get_grid_point_coordinates().transpose((1, 2, 3, 0))
            c_v = gd.cell_cv.sum(0) / 2
            r_gv -= c_v
            lap = Laplace(gd, n=n)
            grad_v = [Gradient(gd, v, n=n) for v in range(3)]
            assert lap.npoints == D * 2 * n + 1
            for m in range(0, 2 * n + 1):
                for ix in range(m + 1):
                    for iy in range(m - ix + 1):
                        iz = m - ix - iy
                        a_g = (r_gv**(ix, iy, iz)).prod(3)
                        if ix + iy + iz == 2 and max(ix, iy, iz) == 2:
                            r = 2.0
                        else:
                            r = 0.0
                        lap.apply(a_g, b_g)
                        e = b_g[n + 1, n + 1, n + 1] - r
                        assert abs(e) < 2e-12, e
                        for v in range(3):
                            grad_v[v].apply(a_g, b_g)
                            if m == 1 and [ix, iy, iz][v] == 1:
                                r = 1
                            else:
                                r = 0
                            e = b_g[n + 1, n + 1, n + 1] - r
                            assert abs(e) < 4e-12, (n, ix, iy, iz, r, v, e)
