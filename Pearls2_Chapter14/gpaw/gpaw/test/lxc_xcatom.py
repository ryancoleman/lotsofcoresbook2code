from __future__ import print_function
import numpy as np
import numpy.random as ra
from gpaw.setup import create_setup
from gpaw.xc import XC
from gpaw.test import equal, gen

if 1:
    for functional in [
        'LDA_X', 'LDA_X+LDA_C_PW', 'LDA_X+LDA_C_VWN', 'LDA_X+LDA_C_PZ',
        'GGA_X_PBE+GGA_C_PBE', 'GGA_X_PBE_R+GGA_C_PBE',
        'GGA_X_B88+GGA_C_P86', 'GGA_X_B88+GGA_C_LYP',
        'GGA_X_FT97_A+GGA_C_LYP'
        ]:
        gen('N', xcname=functional)
        
tolerance = 0.000005 # libxc must reproduce old gpaw energies
# zero Kelvin: in Hartree

reference = { # version 0.9.1
    'LDA_X+LDA_C_PW': 2.28836113207, # 'LDA'
    'GGA_X_PBE+GGA_C_PBE': 2.3366049993, # 'PBE'
    'GGA_X_PBE_R+GGA_C_PBE': 2.34496288319, # 'revPBE'
    }

tolerance_libxc = 0.000001 # libxc must reproduce reference libxc energies

reference_libxc = { # svnversion 5252
    'LDA_X': 1.95030600807,
    'LDA_X+LDA_C_PW': 2.23194461135,
    'LDA_X+LDA_C_VWN': 2.23297429824,
    'LDA_X+LDA_C_PZ': 2.23146045547,
    'GGA_X_PBE+GGA_C_PBE': 2.28208665019,
    'GGA_X_PBE_R+GGA_C_PBE': 2.29201920843,
    'GGA_X_B88+GGA_C_P86': 2.30508027546,
    'GGA_X_B88+GGA_C_LYP': 2.28183010548,
    'GGA_X_FT97_A+GGA_C_LYP': 2.26846048873,
    }

libxc_set = [
    'LDA_X', 'LDA_X+LDA_C_PW', 'LDA_X+LDA_C_VWN', 'LDA_X+LDA_C_PZ',
    'GGA_X_PBE+GGA_C_PBE', 'GGA_X_PBE_R+GGA_C_PBE',
    'GGA_X_B88+GGA_C_P86', 'GGA_X_B88+GGA_C_LYP',
    'GGA_X_FT97_A+GGA_C_LYP'
    ]

x = 0.000001
for xcname in libxc_set:
    ra.seed(8)
    xc = XC(xcname)
    s = create_setup('N', xc)
    ni = s.ni
    nii = ni * (ni + 1) // 2
    D_p = 0.1 * ra.random(nii) + 0.4
    H_p = np.zeros(nii)

    E1 = xc.calculate_paw_correction(s, D_p.reshape(1, -1), H_p.reshape(1, -1))
    dD_p = x * ra.random(nii)
    D_p += dD_p
    dE = np.dot(H_p, dD_p) / x
    E2 = xc.calculate_paw_correction(s, D_p.reshape(1, -1))
    print(xcname, dE, (E2 - E1) / x)
    equal(dE, (E2 - E1) / x, 0.003)

    E2s = xc.calculate_paw_correction(s,
        np.array([0.5 * D_p, 0.5 * D_p]), np.array([H_p, H_p]))
    print(E2, E2s)
    equal(E2, E2s, 1.0e-12)

    if xcname in reference: # compare with old gpaw
        print('A:', E2, reference[xcname])
        equal(E2, reference[xcname], tolerance)

    if xc in reference_libxc: # compare with reference libxc
        print('B:', E2, reference_libxc[xcname])
        equal(E2, reference_libxc[xcname], tolerance)

    D_sp = 0.1 * ra.random((2, nii)) + 0.2
    H_sp = np.zeros((2, nii))

    E1 = xc.calculate_paw_correction(s, D_sp, H_sp)
    dD_sp = x * ra.random((2, nii))
    D_sp += dD_sp
    dE = np.dot(H_sp.ravel(), dD_sp.ravel()) / x
    E2 = xc.calculate_paw_correction(s, D_sp, H_sp)
    print(dE, (E2 - E1) / x)
    equal(dE, (E2 - E1) / x, 0.005)
