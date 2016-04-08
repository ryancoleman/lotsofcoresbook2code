# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

from math import pi, sqrt

import numpy as np

from gpaw.utilities.blas import axpy, gemm, gemv, gemmdot
from gpaw import extra_parameters
from gpaw.gaunt import gaunt
from gpaw.spherical_harmonics import nablarlYL

# load points and weights for the angular integration
from gpaw.sphere.lebedev import Y_nL, R_nv, weight_n


"""
                           3
             __   dn       __   __    dY
   __  2    \       L  2  \    \        L  2
  (\/n) = (  ) Y  --- )  + ) (  )  n  --- )
            /__ L dr      /__  /__  L dr
                                        v
             L            v=1    L


        dY
          L
  A   = --- r
   Lv   dr
          v

"""
# A_nvL is defined as above, n is an expansion point index (50 Lebedev points).
rnablaY_nLv = np.empty((len(R_nv), 25, 3))
for rnablaY_Lv, Y_L, R_v in zip(rnablaY_nLv, Y_nL, R_nv):
    for l in range(5):
        for L in range(l**2, (l + 1)**2):
            rnablaY_Lv[L] = nablarlYL(L, R_v)  - l * R_v * Y_L[L]


class PAWXCCorrection:
    def __init__(self,
                 w_jg,  # all-lectron partial waves
                 wt_jg, # pseudo partial waves
                 nc_g,  # core density
                 nct_g, # smooth core density
                 rgd,   # radial grid descriptor
                 jl,    # ?
                 lmax,  # maximal angular momentum to consider
                 Exc0,  # xc energy of reference atom
                 phicorehole_g, # ?
                 fcorehole,     # ?
                 tauc_g,   # kinetic core energy array
                 tauct_g   # pseudo kinetic core energy array
                 ):
        self.Exc0 = Exc0
        self.Lmax = (lmax + 1)**2
        self.rgd = rgd
        self.dv_g = rgd.dv_g
        #self.nspins = nspins
        self.Y_nL = Y_nL[:, :self.Lmax]
        self.rnablaY_nLv = rnablaY_nLv[:, :self.Lmax]
        self.ng = ng = len(nc_g)
        self.phi_jg = w_jg
        self.phit_jg = wt_jg
        
        self.jlL = [(j, l, l**2 + m) for j, l in jl for m in range(2 * l + 1)]
        self.ni = ni = len(self.jlL)
        self.nj = nj = len(jl)
        self.nii = nii = ni * (ni + 1) // 2
        njj = nj * (nj + 1) // 2

        self.tauc_g = tauc_g
        self.tauct_g = tauct_g
        self.tau_npg = None
        self.taut_npg = None

        B_Lqp = np.zeros((self.Lmax, njj, nii))
        p = 0
        i1 = 0
        for j1, l1, L1 in self.jlL:
            for j2, l2, L2 in self.jlL[i1:]:
                if j1 < j2:
                    q = j2 + j1 * nj - j1 * (j1 + 1) // 2
                else:
                    q = j1 + j2 * nj - j2 * (j2 + 1) // 2
                B_Lqp[:, q, p] = gaunt[L1, L2, :self.Lmax]
                p += 1
            i1 += 1
        self.B_pqL = B_Lqp.T.copy()

        #
        self.n_qg = np.zeros((njj, ng))
        self.nt_qg = np.zeros((njj, ng))
        q = 0
        for j1, l1 in jl:
            for j2, l2 in jl[j1:]:
                #rl1l2 = rgd.r_g**(l1 + l2)
                self.n_qg[q] = w_jg[j1] * w_jg[j2]
                self.nt_qg[q] = wt_jg[j1] * wt_jg[j2]
                q += 1

        self.nc_g = nc_g
        self.nct_g = nct_g

        if fcorehole != 0.0:
            self.nc_corehole_g = fcorehole * phicorehole_g**2 / (4 * pi)
        else:
            self.nc_corehole_g = None

    def four_phi_integrals(self, D_sp, fxc):
        """Calculate four-phi integrals.

        The density is given by the density matrix ``D_sp`` in packed(pack)
        form, and the resulting rank-four tensor is also returned in
        packed format. ``fxc`` is a radial object???
        """

        ns, nii = D_sp.shape

        assert ns == 1# and not self.xc.get_functional().gga

        D_p = D_sp[0]
        D_Lq = np.dot(self.B_pqL.T, D_p)

        # Expand all-electron density in spherical harmonics:
        n_qg = self.n_qg
        n_Lg = np.dot(D_Lq, n_qg)
        n_Lg[0] += self.nc_g * sqrt(4 * pi)

        # Expand pseudo electron density in spherical harmonics:
        nt_qg = self.nt_qg
        nt_Lg = np.dot(D_Lq, nt_qg)
        nt_Lg[0] += self.nct_g * sqrt(4 * pi)

        # Allocate array for result:
        J_pp = np.zeros((nii, nii))

        # Loop over 50 points on the sphere surface:
        for w, Y_L in zip(weight_n, self.Y_nL):
            B_pq = np.dot(self.B_pqL, Y_L)

            fxcdv = fxc(np.dot(Y_L, n_Lg)) * self.dv_g
            dn2_qq = np.inner(n_qg * fxcdv, n_qg)

            fxctdv = fxc(np.dot(Y_L, nt_Lg)) * self.dv_g
            dn2_qq -= np.inner(nt_qg * fxctdv, nt_qg)

            J_pp += w * np.dot(B_pq, np.inner(dn2_qq, B_pq))

        return J_pp
