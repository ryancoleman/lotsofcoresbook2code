import numpy as np

# Test that numpy.dot behaves as expected, i.e.
# [A . B]_ijpq = sum_k  A_ijk * B_pkq
# Product sum is over last dimension of A and second-to-last dimension of B
#
# See numpy.tensordot for ultimate flexibility in choosing the pruduct-sum axes

# make "random" input arrays
A = np.arange(6 * 2 * 5).reshape(6, 2, 5) - 3.
B = np.arange(3 * 5 * 4).reshape(3, 5, 4) + 5.

# built-in dot product
AB1 = np.dot(A, B)

# manual dot product
AB2 = np.empty(A.shape[:-1] + B.shape[:-2] + (B.shape[-1],))
for i in range(AB2.shape[0]):
  for j in range(AB2.shape[1]):
     for p in range(AB2.shape[2]):
       for q in range(AB2.shape[3]):
         AB2[i, j, p, q] = np.sum(A[i, j, :] * B[p, :, q])

# test for consistency
#print AB1 == AB2
assert np.all(AB1 == AB2)
