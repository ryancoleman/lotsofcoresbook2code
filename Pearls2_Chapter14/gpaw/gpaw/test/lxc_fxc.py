from __future__ import print_function
from gpaw.xc.libxc import LibXC
from math import pi
import numpy as np

nspins = 1
for name in [
    'LDA', 'PBE', 'revPBE', 'RPBE',
    'LDA_X', 'GGA_X_PBE_R', 'GGA_X_RPBE',
    'LDA_C_PW',
    ]:
    xc = LibXC(name)
    xc.initialize(nspins)
    libxc = xc.xc
    lxc_fxc = libxc.calculate_fxc_spinpaired
    lxc_fxc_fd = libxc.calculate_fxc_fd_spinpaired
    na = 2.0
    if nspins == 2:
        nb = 1.0
    else:
        nb = 0.0
    print(na, nb)
    if (nb > 0.0): assert (nspins == 2)
    if nspins == 2:
        sigma0 = 2.0 # (0.0, 1.0, 1.0)
        sigma1 = 2.0
        sigma2 = 5.0 # (1.0, 2.0, 0.0)
    else:
        sigma0 = 2.0 # (0.0, 1.0, 1.0)
        sigma1 = 0.0
        sigma2 = 0.0
    taua=(3.*pi**2)**(2./3.)*na**(5./3.)/2.*sigma0
    taub=(3.*pi**2)**(2./3.)*nb**(5./3.)/2.*sigma2
    if ((sigma1 > 0.0) or (sigma2 > 0.0)): assert (nspins == 2)

    na = np.array([na])
    sigma0 = np.array([sigma0])

    dvdn = np.zeros((1))
    dvdnda2 = np.zeros((6))
    dvda2da2 = np.zeros_like(dvdnda2)

    dvdn_N = np.zeros_like(dvdn)
    dvdnda2_N = np.zeros_like(dvdnda2)
    dvda2da2_N = np.zeros_like(dvdnda2)

    lxc_fxc(na, dvdn, sigma0, dvdnda2, dvda2da2)
    lxc_fxc_fd(na, dvdn_N, sigma0, dvdnda2_N, dvda2da2_N)

    error = [0.0, 'exact']
    for E in [
        ('dvdn', dvdn[0], dvdn_N[0]),
        ('dvdnda2', dvdnda2[0], dvdnda2_N[0]),
        ('dvda2da2', dvda2da2[0], dvda2da2_N[0]),
        ]:
        for e in E[2:]:
            de = abs(e - E[1])
            if de > error[0]:
                error[0] = de
                error[1] = E[0]
    print(name, error[0], error[1])
    assert error[0] < 5.0e-7
