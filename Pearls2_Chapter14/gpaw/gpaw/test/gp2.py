from __future__ import print_function
import numpy as np
from gpaw.grid_descriptor import GridDescriptor
from gpaw.localized_functions import create_localized_functions
from gpaw.spline import Spline

s=Spline(0, 1.2, [1, 0.6, 0.1, 0.0])
a = 4.0
n = 24
gd = GridDescriptor((n, n, n), (a, a, a))
print(gd.get_boxes((0, 0, 0), 1.2, 0))
if 0:
    p = create_localized_functions([s], gd, (0.0, 0.0, 0.0), cut=True)
    a = np.zeros((n, n, n))
    p.add(a, np.array([2.0]))
    print(a[1,0])
