# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

from math import sqrt, cos, sin

import numpy as np

from gpaw.spherical_harmonics import Y



s = sqrt(0.5)
t = sqrt(3) / 3
# Points on the unit sphere:
sphere_lm = [ \
    np.array([(1, 0, 0)]), # s
    np.array([(1, 0, 0), (0, 1, 0), (0, 0, 1)]), # p
    np.array([(s, s, 0), (0, s, s), (s, 0, s), (1, 0, 0), (0, 0, 1)]), # d
    np.array([(s, s, 0), (0, s, s), (s, 0, s),
              (1, 0, 0), (s, -s, 0),
              (t, -t, t), (t, t, -t)])] # f

def Y_matrix(l, U_vv):
    """YMatrix(l, symmetry) -> matrix.

    For 2*l+1 points on the unit sphere (m1=0,...,2*l) calculate the
    value of Y_lm2 for m2=0,...,2*l.  The points are those from the
    list sphere_lm[l] rotated as described by U_vv."""
    
    Y_mm = np.zeros((2 * l + 1, 2 * l + 1))
    for m1, point in enumerate(sphere_lm[l]):
        x, y, z = np.dot(point, U_vv)
        for m2 in range(2 * l + 1):
            L = l**2 + m2
            Y_mm[m1, m2] = Y(L, x, y, z)
    return Y_mm

identity = np.eye(3)
iY_lmm = [np.linalg.inv(Y_matrix(l, identity)) for l in range(4)]
         

def rotation(l, U_vv):
    """Rotation(l, symmetry) -> transformation matrix.

    Find the transformation from Y_lm1 to Y_lm2."""
    
    return np.dot(iY_lmm[l], Y_matrix(l, U_vv))


def Y_rotation(l, angle):
    Y_mm = np.zeros((2 * l + 1, 2 * l + 1))
    sn = sin(angle)
    cs = cos(angle)
    for m1, point in enumerate(sphere_lm[l]):
        x = point[0]
        y = cs * point[1] - sn * point[2]
        z = sn * point[1] + cs * point[2]
        for m2 in range(2 * l + 1):
            L = l**2 + m2
            Y_mm[m1, m2] = Y(L, x, y, z)
    return Y_mm


del s, identity
