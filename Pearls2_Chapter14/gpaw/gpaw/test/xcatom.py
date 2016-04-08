from __future__ import print_function
import numpy as np
import numpy.random as ra
from gpaw.setup import create_setup
from gpaw.xc import XC
from gpaw.test import equal


x = 0.000001
ra.seed(8)
for xc in ['LDA', 'PBE']:
    print(xc)
    xc = XC(xc)
    s = create_setup('N', xc)
    ni = s.ni
    nii = ni * (ni + 1) // 2
    D_p = 0.1 * ra.random(nii) + 0.2
    H_p = np.zeros(nii)

    E = xc.calculate_paw_correction(s, D_p.reshape(1, -1), H_p.reshape(1, -1))
    dD_p = x * ra.random(nii)
    dE = np.dot(H_p, dD_p) / x
    D_p += dD_p
    Ep = xc.calculate_paw_correction(s, D_p.reshape(1, -1))
    D_p -= 2 * dD_p
    Em = xc.calculate_paw_correction(s, D_p.reshape(1, -1))
    print(dE, dE - 0.5 * (Ep - Em) / x)
    equal(dE, 0.5 * (Ep - Em) / x, 1e-6)

    Ems = xc.calculate_paw_correction(s, np.array([0.5 * D_p, 0.5 * D_p]))
    print(Em - Ems)
    equal(Em, Ems, 1.0e-12)

    D_sp = 0.1 * ra.random((2, nii)) + 0.2
    H_sp = np.zeros((2, nii))

    E = xc.calculate_paw_correction(s, D_sp, H_sp)
    dD_sp = x * ra.random((2, nii))
    dE = np.dot(H_sp.ravel(), dD_sp.ravel()) / x
    D_sp += dD_sp
    Ep = xc.calculate_paw_correction(s, D_sp)
    D_sp -= 2 * dD_sp
    Em = xc.calculate_paw_correction(s, D_sp)
    print(dE, dE - 0.5 * (Ep - Em) / x)
    equal(dE, 0.5 * (Ep - Em) / x, 1e-6)
