from __future__ import print_function
import numpy as np

pos = np.asarray([1.450 * a for a in range(6)])
n = len(pos)
d2_ii = np.zeros((n, n))
for i1, i2 in np.ndindex(n, n):
    d2_ii[i1, i2] = ((pos[i1]-pos[i2])**2).sum()

d_ii = np.sqrt(d2_ii)
U = 11.26
t_m = -2.4
t_lm = -1.0
# Two particle interaction
V = U / np.sqrt(1.0 + 0.6117*d2_ii)
V[0,-1] = V[-1,0] = 0.0
V[0,0] = V[-1,-1] = 0.0
# One particle part
mask_c = d_ii < 1.5 
mask_c *= d_ii > 0.0 
mask_c = mask_c.astype(int)
nbf = len(d_ii)
h = t_m * np.ones((nbf, nbf)) * mask_c
h[0, 1] = h[1, 0] = t_lm
h[-2,-1] = h[-1,-2] = t_lm

# H_lead
h1 = np.zeros((2,2), complex)
t_l = -20.0
h1[0,1] = t_l
h1[1,0] = t_l
nbf = len(h)

# H_scat
H = np.zeros((nbf+2, nbf+2), complex) 
H[0,1] = H[1,0] = t_l
H[-2,-1] = H[-1,-2] = t_l
H[1:-1, 1:-1] = h

# Hartree potential of the ions (Z=1)
ion_shift = np.zeros((n, n))
for i in range(n):
    ion_shift[i, i] += -0.5 * V[i, i]
    for k in range(n):
        if k!=i:
            ion_shift[i, i] += -1.0 * V[i, k]

H = H.astype(complex)
V = V.astype(complex)
ion_shift = ion_shift.astype(complex)

if __name__=='__main__':
    print(np.around(H.real, 2))
    print(np.around(V.real ,2))
    print(np.around(ion_shift.diagonal().real, 2))

