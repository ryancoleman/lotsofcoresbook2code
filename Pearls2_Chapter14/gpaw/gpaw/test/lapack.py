from __future__ import print_function
import numpy as np

from gpaw.test import equal
from gpaw.utilities.lapack import sqrt_matrix

# check sqrt of a matrix

A = [[20, 4], [4, 1]]
a = [[4.4, .8], [.8, .6]]
A = np.array(A, float)
print('A=', A)
a = np.array(a)
b = sqrt_matrix(A)
print('sqrt(A)=', b)
equal(((a-b)**2).sum(), 0, 1.e-12)
