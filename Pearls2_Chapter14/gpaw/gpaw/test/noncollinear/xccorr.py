from __future__ import print_function
import numpy as np
import numpy.random as ra
from gpaw.setup import create_setup
from gpaw.xc.noncollinear import NonCollinearFunctional
from gpaw.xc import XC
from gpaw.test import equal
from gpaw.utilities import pack


x = 0.0000001
ra.seed(8)
for xc in ['LDA']:#, 'PBE']:
    print(xc)
    xc = XC(xc)
    s = create_setup('N', xc)
    ni = s.ni
    D_sp = np.array([pack(np.outer(P_i, P_i))
                     for P_i in 0.2 * ra.random((2, ni)) - 0.1])
    D_sp[1] += D_sp[0]

    nii = ni * (ni + 1) // 2

    E = xc.calculate_paw_correction(s, D_sp)

    xc = NonCollinearFunctional(xc)

    Dnc_sp = np.zeros((4, nii))
    Dnc_sp[0] = D_sp.sum(0)
    Dnc_sp[3] = D_sp[0] - D_sp[1]
    Enc = xc.calculate_paw_correction(s, Dnc_sp)
    print(E, E-Enc)
    assert abs(E - Enc) < 1e-11

    Dnc_sp[1] = 0.1 * Dnc_sp[3]
    Dnc_sp[2] = 0.2 * Dnc_sp[3]
    Dnc_sp[3] *= (1 - 0.1**2 - 0.2**2)**0.5
    H_sp = 0 * Dnc_sp
    Enc = xc.calculate_paw_correction(s, Dnc_sp, H_sp)
    print(E, E-Enc)
    assert abs(E - Enc) < 1e-11

    dD_sp = x * ra.random((4, nii))
    dE = np.dot(H_sp.ravel(), dD_sp.ravel()) / x
    Dnc_sp += dD_sp
    Ep = xc.calculate_paw_correction(s, Dnc_sp)
    Dnc_sp -= 2 * dD_sp
    Em = xc.calculate_paw_correction(s, Dnc_sp)
    print(dE, dE - 0.5 * (Ep - Em) / x)
    assert abs(dE - 0.5 * (Ep - Em) / x) < 2e-8
        
