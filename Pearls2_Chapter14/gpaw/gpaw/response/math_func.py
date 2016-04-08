from math import sqrt, pi

import numpy as np
from scipy.special.specfun import sphj

from gpaw.utilities.blas import gemmdot
from gpaw.gaunt import gaunt as G_LLL
from gpaw.spherical_harmonics import Y


def delta_function(x0, dx, Nx, sigma):

    deltax = np.zeros(Nx)
    for i in range(Nx):
        deltax[i] = np.exp(-(i * dx - x0)**2/(4. * sigma))
    return deltax / (2. * sqrt(pi * sigma))


def hilbert_transform(specfunc_wGG, w_w, Nw, dw, eta, fullresponse=False):

    NwS = specfunc_wGG.shape[0]
    tmp_ww = np.zeros((Nw, NwS), dtype=complex)
    ww_w = np.linspace(0., (NwS-1)*dw, NwS)

    for iw in range(Nw):
        if fullresponse is False:
            tmp_ww[iw] = 1. / (w_w[iw] - ww_w + 1j*eta) - 1. / (w_w[iw] + ww_w + 1j*eta)
        else:
            tmp_ww[iw] = 1. / (w_w[iw] - ww_w + 1j*eta) - 1. / (w_w[iw] + ww_w - 1j*eta)

    chi0_wGG = gemmdot(tmp_ww, specfunc_wGG, beta = 0.)

    return chi0_wGG * dw


def two_phi_planewave_integrals(k_Gv, setup=None, Gstart=0, Gend=None,
                                rgd=None, phi_jg=None,
                                phit_jg=None,l_j=None):
    """Calculate PAW-correction matrix elements with planewaves.

    ::
    
      /  _       _   ik.r     _     ~   _   ik.r ~   _
      | dr [phi (r) e    phi (r) - phi (r) e    phi (r)]
      /        1            2         1            2

                        ll    -  /     2                      ~       ~
      = 4 * pi \sum_lm  i  Y (k) | dr r  [ phi (r) phi (r) - phi (r) phi (r) j (kr)
                            lm   /            1       2         1       2     ll

           /
        * | d\Omega Y     Y     Y
          /          l1m1  l2m2  lm
          
    """

    if Gend is None:
        Gend = len(k_Gv)
        
    if setup is not None:
        rgd = setup.rgd
        l_j = setup.l_j
        # Obtain the phi_j and phit_j
        phi_jg = []
        phit_jg = []
        rcut2 = 2 * max(setup.rcut_j)
        gcut2 = rgd.ceil(rcut2)
        for phi_g, phit_g in zip(setup.data.phi_jg, setup.data.phit_jg):
            phi_g = phi_g.copy()
            phit_g = phit_g.copy()
            phi_g[gcut2:] = phit_g[gcut2:] = 0.
            phi_jg.append(phi_g)
            phit_jg.append(phit_g)
    else:
        assert rgd is not None
        assert phi_jg is not None
        assert l_j is not None

    ng = rgd.N
    r_g = rgd.r_g
    dr_g = rgd.dr_g

    # Construct L (l**2 + m) and j (nl) index
    L_i = []
    j_i = []
    for j, l in enumerate(l_j):
        for m in range(2 * l + 1):
            L_i.append(l**2 + m)
            j_i.append(j)
    ni = len(L_i)
    nj = len(l_j)
    lmax = max(l_j) * 2 + 1

    if setup is not None:
        assert ni == setup.ni and nj == setup.nj

    # Initialize
    npw = k_Gv.shape[0]
    R_jj = np.zeros((nj, nj))
    R_ii = np.zeros((ni, ni))
    phi_Gii = np.zeros((npw, ni, ni), dtype=complex)
    j_lg = np.zeros((lmax, ng))

    # Store (phi_j1 * phi_j2 - phit_j1 * phit_j2 ) for further use
    tmp_jjg = np.zeros((nj, nj, ng))
    for j1 in range(nj):
        for j2 in range(nj):
            tmp_jjg[j1, j2] = (phi_jg[j1] * phi_jg[j2] -
                               phit_jg[j1] * phit_jg[j2])

    # Loop over G vectors
    for iG in range(Gstart, Gend):
        kk = k_Gv[iG] 
        k = np.sqrt(np.dot(kk, kk)) # calculate length of q+G

        # Calculating spherical bessel function
        for g, r in enumerate(r_g):
            j_lg[:, g] = sphj(lmax - 1,  k * r)[1]

        for li in range(lmax):
            # Radial part 
            for j1 in range(nj):
                for j2 in range(nj): 
                    R_jj[j1, j2] = np.dot(r_g**2*dr_g,
                                          tmp_jjg[j1, j2] * j_lg[li])

            for mi in range(2 * li + 1):
                # Angular part
                for i1 in range(ni):
                    L1 = L_i[i1]
                    j1 = j_i[i1]
                    for i2 in range(ni):
                        L2 = L_i[i2]
                        j2 = j_i[i2]
                        R_ii[i1, i2] = G_LLL[L1, L2, li**2+mi]  * R_jj[j1, j2]

                phi_Gii[iG] += R_ii * Y(li**2 + mi,
                                        kk[0]/k, kk[1]/k, kk[2]/k) * (-1j)**li
    
    phi_Gii *= 4 * pi

    return phi_Gii.reshape(npw, ni*ni)
