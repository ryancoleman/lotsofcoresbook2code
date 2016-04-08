# -*- coding: utf-8
# In an attempt to appease epydoc and still have readable docstrings,
# the vertical bar | is represented by u'\u2758' in this module.
"""This module defines Coulomb and XC kernels for the response model.
"""

import numpy as np
#from math import sqrt, pi, sin, cos, exp
from gpaw.utilities.blas import gemmdot
from gpaw.xc import XC
from gpaw.sphere.lebedev import weight_n, R_nv
from gpaw.mpi import world, rank, size
from ase.dft.kpoints import monkhorst_pack


def truncated_coulomb(pd):
    """ Simple truncation of Coulomb kernel along z
    Rozzi, C., Varsano, D., Marini, A., Gross, E., & Rubio, A. (2006).
    Exact Coulomb cutoff technique for supercell calculations.
    Physical Review B, 73(20), 205119. doi:10.1103/PhysRevB.73.205119
    """

    qG = pd.get_reciprocal_vectors(add_q=True)
    if pd.kd.gamma:  # Set small finite q to handle divergence
            qG += [0.0000001, 0, 0]
    nG = len(qG)
    L = pd.gd.cell_cv[2, 2]
    R = L / 2.  # Truncation length is half of unit cell
    qG_par = ((qG[:, 0])**2 + (qG[:, 1]**2))**0.5
    qG_z = qG[:, 2]
    
    K_G = 4 * np.pi / (qG**2).sum(axis=1)

    # K_G *= 1. + np.exp(-qG_par * R) * (qG_z / qG_par * np.sin(qG_z * R)\
    #                                     - np.cos(qG_z * R))
    # sin(qG_z * R) = 0 when R = L/2

    K_G *= 1. - np.exp(-qG_par * R) * np.cos(qG_z * R)
    
    K_G **= 0.5

    return K_G.astype(complex)


def calculate_Kxc(pd, nt_sG, R_av, setups, D_asp, functional='ALDA',
                  density_cut=None):
    """ALDA kernel"""

    gd = pd.gd
    npw = pd.ngmax
    nG = pd.gd.N_c
    vol = pd.gd.volume
    bcell_cv = np.linalg.inv(pd.gd.cell_cv)
    G_Gv = pd.get_reciprocal_vectors()
    
    # The soft part
    #assert np.abs(nt_sG[0].shape - nG).sum() == 0
    if functional == 'ALDA_X':
        x_only = True
        A_x = -3. / 4. * (3. / np.pi)**(1. / 3.)
        nspins = len(nt_sG)
        assert nspins in [1, 2]
        fxc_sg = nspins**(1. / 3.) * 4. / 9. * A_x * nt_sG**(-2. / 3.)
    else:
        assert len(nt_sG) == 1
        x_only = False
        fxc_sg = np.zeros_like(nt_sG)
        xc = XC(functional[1:])
        xc.calculate_fxc(gd, nt_sG, fxc_sg)

    if density_cut is not None:
        fxc_sg[np.where(nt_sG * len(nt_sG) < density_cut)] = 0.0

    # FFT fxc(r)
    nG0 = nG[0] * nG[1] * nG[2]
    tmp_sg = [np.fft.fftn(fxc_sg[s]) * vol / nG0 for s in range(len(nt_sG))]

    r_vg = gd.get_grid_point_coordinates()
    Kxc_sGG = np.zeros((len(fxc_sg), npw, npw), dtype=complex)
    for s in range(len(fxc_sg)):
        for iG, iQ in enumerate(pd.Q_qG[0]):
            iQ_c = (np.unravel_index(iQ, nG) + nG // 2) % nG - nG // 2
            for jG, jQ in enumerate(pd.Q_qG[0]):
                jQ_c = (np.unravel_index(jQ, nG) + nG // 2) % nG - nG // 2
                ijQ_c = (iQ_c - jQ_c)
                if (abs(ijQ_c) < nG // 2).all():
                    Kxc_sGG[s, iG, jG] = tmp_sg[s][tuple(ijQ_c)]

    # The PAW part
    KxcPAW_sGG = np.zeros_like(Kxc_sGG)
    dG_GGv = np.zeros((npw, npw, 3))
    for v in range(3):
        dG_GGv[:, :, v] = np.subtract.outer(G_Gv[:, v], G_Gv[:, v])

    for a, setup in enumerate(setups):
        if rank == a % size:
            rgd = setup.xc_correction.rgd
            n_qg = setup.xc_correction.n_qg
            nt_qg = setup.xc_correction.nt_qg
            nc_g = setup.xc_correction.nc_g
            nct_g = setup.xc_correction.nct_g
            Y_nL = setup.xc_correction.Y_nL
            dv_g = rgd.dv_g

            D_sp = D_asp[a]
            B_pqL = setup.xc_correction.B_pqL
            D_sLq = np.inner(D_sp, B_pqL.T)
            nspins = len(D_sp)

            f_sg = rgd.empty(nspins)
            ft_sg = rgd.empty(nspins)

            n_sLg = np.dot(D_sLq, n_qg)
            nt_sLg = np.dot(D_sLq, nt_qg)

            # Add core density
            n_sLg[:, 0] += np.sqrt(4. * np.pi) / nspins * nc_g
            nt_sLg[:, 0] += np.sqrt(4. * np.pi) / nspins * nct_g

            coefatoms_GG = np.exp(-1j * np.inner(dG_GGv, R_av[a]))
            for n, Y_L in enumerate(Y_nL):
                w = weight_n[n]
                f_sg[:] = 0.0
                n_sg = np.dot(Y_L, n_sLg)
                if x_only:
                    f_sg = nspins * (4 / 9.) * A_x * (nspins * n_sg)**(-2 / 3.)
                else:
                    xc.calculate_fxc(rgd, n_sg, f_sg)

                ft_sg[:] = 0.0
                nt_sg = np.dot(Y_L, nt_sLg)
                if x_only:
                    ft_sg = nspins * (4 / 9.) * (A_x
                                                 * (nspins * nt_sg)**(-2 / 3.))
                else:
                    xc.calculate_fxc(rgd, nt_sg, ft_sg)
                for i in range(len(rgd.r_g)):
                    coef_GG = np.exp(-1j * np.inner(dG_GGv, R_nv[n])
                                     * rgd.r_g[i])
                    for s in range(len(f_sg)):
                        KxcPAW_sGG[s] += w * np.dot(coef_GG,
                                                    (f_sg[s, i] -
                                                     ft_sg[s, i])
                                                    * dv_g[i]) * coefatoms_GG

    world.sum(KxcPAW_sGG)
    Kxc_sGG += KxcPAW_sGG

    return Kxc_sGG / vol
