from math import sqrt, pi

import numpy as np

from gpaw.xc.gga import GGA
from gpaw.utilities.blas import axpy
from gpaw.sphere.lebedev import weight_n


class MGGA(GGA):
    orbital_dependent = True

    def __init__(self, kernel):
        """Meta GGA functional."""
        GGA.__init__(self, kernel)

    def set_grid_descriptor(self, gd):
        GGA.set_grid_descriptor(self, gd)

    def get_setup_name(self):
        return 'PBE'

    def initialize(self, density, hamiltonian, wfs, occupations):
        self.wfs = wfs
        self.tauct = density.get_pseudo_core_kinetic_energy_density_lfc()
        self.tauct_G = None
        self.dedtaut_sG = None
        self.restrict = hamiltonian.restrict
        self.interpolate = density.interpolate

    def set_positions(self, spos_ac):
        self.tauct.set_positions(spos_ac)
        if self.tauct_G is None:
            self.tauct_G = self.wfs.gd.empty()
        self.tauct_G[:] = 0.0
        self.tauct.add(self.tauct_G)

    def calculate_gga(self, e_g, nt_sg, v_sg, sigma_xg, dedsigma_xg):
        try:
            taut_sG = self.wfs.calculate_kinetic_energy_density()
        except RuntimeError:
            nspins = self.wfs.nspins
            # Initialize with von Weizsaecker kinetic energy density
            taut_sG = self.wfs.gd.empty((nspins))
            gradn_g = self.gd.empty()
            for s in range(nspins):
                taut_g = self.gd.zeros()
                for v in range(3):
                    self.grad_v[v](nt_sg[s], gradn_g)
                    axpy(0.125, gradn_g**2, taut_g)
                ntinv_g = 0. * taut_g
                nt_ok = np.where(nt_sg[s] > 1e-7)
                ntinv_g[nt_ok] = 1.0 / nt_sg[s][nt_ok]
                taut_g *= ntinv_g
                self.restrict(taut_g, taut_sG[s])

        taut_sg = np.empty_like(nt_sg)
        for taut_G, taut_g in zip(taut_sG, taut_sg):
            taut_G += 1.0 / self.wfs.nspins * self.tauct_G
            self.interpolate(taut_G, taut_g)

        dedtaut_sg = np.empty_like(nt_sg)

        self.kernel.calculate(e_g, nt_sg, v_sg, sigma_xg, dedsigma_xg,
                              taut_sg, dedtaut_sg)
        self.dedtaut_sG = self.wfs.gd.empty(self.wfs.nspins)
        self.ekin = 0.0
        for s in range(self.wfs.nspins):
            self.restrict(dedtaut_sg[s], self.dedtaut_sG[s])
            self.ekin -= self.wfs.gd.integrate(
                self.dedtaut_sG[s] * (taut_sG[s] -
                                      self.tauct_G / self.wfs.nspins))

    def apply_orbital_dependent_hamiltonian(self, kpt, psit_xG,
                                            Htpsit_xG, dH_asp):
        self.wfs.apply_mgga_orbital_dependent_hamiltonian(
            kpt, psit_xG,
            Htpsit_xG, dH_asp,
            self.dedtaut_sG[kpt.s])

    def calculate_paw_correction(self, setup, D_sp, dEdD_sp=None,
                                 addcoredensity=True, a=None):
        assert not hasattr(self, 'D_sp')
        self.D_sp = D_sp
        self.n = 0
        self.ae = True
        self.c = setup.xc_correction
        self.dEdD_sp = dEdD_sp

        if self.c.tau_npg is None:
            self.c.tau_npg, self.c.taut_npg = self.initialize_kinetic(self.c)

        E = GGA.calculate_paw_correction(self, setup, D_sp, dEdD_sp,
                                         addcoredensity, a)
        del self.D_sp, self.n, self.ae, self.c, self.dEdD_sp
        return E

    def calculate_gga_radial(self, e_g, n_sg, v_sg, sigma_xg, dedsigma_xg):
        nspins = len(n_sg)
        if self.ae:
            tau_pg = self.c.tau_npg[self.n]
            tauc_g = self.c.tauc_g / (sqrt(4 * pi) * nspins)
            sign = 1.0
        else:
            tau_pg = self.c.taut_npg[self.n]
            tauc_g = self.c.tauct_g / (sqrt(4 * pi) * nspins)
            sign = -1.0
        tau_sg = np.dot(self.D_sp, tau_pg) + tauc_g
        dedtau_sg = np.empty_like(tau_sg)
        self.kernel.calculate(e_g, n_sg, v_sg, sigma_xg, dedsigma_xg,
                              tau_sg, dedtau_sg)
        if self.dEdD_sp is not None:
            self.dEdD_sp += (sign * weight_n[self.n] *
                             np.inner(dedtau_sg * self.c.rgd.dv_g, tau_pg))
        self.n += 1
        if self.n == len(weight_n):
            self.n = 0
            self.ae = False

    def calculate_spherical(self, rgd, n_sg, v_sg):
        raise NotImplementedError

    def add_forces(self, F_av):
        dF_av = self.tauct.dict(derivative=True)
        self.tauct.derivative(self.dedtaut_sG.sum(0), dF_av)
        for a, dF_v in dF_av.items():
            F_av[a] += dF_v[0] / self.wfs.nspins

    def estimate_memory(self, mem):
        bytecount = self.wfs.gd.bytecount()
        mem.subnode('MGGA arrays', (1 + self.wfs.nspins) * bytecount)

    def initialize_kinetic(self, xccorr):
        nii = xccorr.nii
        nn = len(xccorr.rnablaY_nLv)
        ng = len(xccorr.phi_jg[0])

        tau_npg = np.zeros((nn, nii, ng))
        taut_npg = np.zeros((nn, nii, ng))
        self.create_kinetic(xccorr, nn, xccorr.phi_jg, tau_npg)
        self.create_kinetic(xccorr, nn, xccorr.phit_jg, taut_npg)
        return tau_npg, taut_npg

    def create_kinetic(self, x, ny, phi_jg, tau_ypg):
        """Short title here.

        kinetic expression is::

                                             __         __
          tau_s = 1/2 Sum_{i1,i2} D(s,i1,i2) \/phi_i1 . \/phi_i2 +tauc_s

        here the orbital dependent part is calculated::

          __         __
          \/phi_i1 . \/phi_i2 =
                      __    __
                      \/YL1.\/YL2 phi_j1 phi_j2 +YL1 YL2 dphi_j1 dphi_j2
                                                         ------  ------
                                                           dr     dr
          __    __
          \/YL1.\/YL2 [y] = Sum_c A[L1,c,y] A[L2,c,y] / r**2

        """
        nj = len(phi_jg)
        dphidr_jg = np.zeros(np.shape(phi_jg))
        for j in range(nj):
            phi_g = phi_jg[j]
            x.rgd.derivative(phi_g, dphidr_jg[j])

        # Second term:
        for y in range(ny):
            i1 = 0
            p = 0
            Y_L = x.Y_nL[y]
            for j1, l1, L1 in x.jlL:
                for j2, l2, L2 in x.jlL[i1:]:
                    c = Y_L[L1] * Y_L[L2]
                    temp = c * dphidr_jg[j1] * dphidr_jg[j2]
                    tau_ypg[y, p, :] += temp
                    p += 1
                i1 += 1
        ##first term
        for y in range(ny):
            i1 = 0
            p = 0
            rnablaY_Lv = x.rnablaY_nLv[y, :x.Lmax]
            Ax_L = rnablaY_Lv[:, 0]
            Ay_L = rnablaY_Lv[:, 1]
            Az_L = rnablaY_Lv[:, 2]
            for j1, l1, L1 in x.jlL:
                for j2, l2, L2 in x.jlL[i1:]:
                    temp = (Ax_L[L1] * Ax_L[L2] + Ay_L[L1] * Ay_L[L2]
                            + Az_L[L1] * Az_L[L2])
                    temp *= phi_jg[j1] * phi_jg[j2]
                    temp[1:] /= x.rgd.r_g[1:]**2
                    temp[0] = temp[1]
                    tau_ypg[y, p, :] += temp
                    p += 1
                i1 += 1
        tau_ypg *= 0.5

        return


class PurePython2DMGGAKernel:
    def __init__(self, name, pars=None):
        self.name = name
        self.pars = pars
        self.type = 'MGGA'
        assert self.pars is not None

    def calculate(self, e_g, n_sg, dedn_sg,
                  sigma_xg, dedsigma_xg,
                  tau_sg, dedtau_sg):

        e_g[:] = 0.
        dedsigma_xg[:] = 0.
        dedtau_sg[:] = 0.

        # spin-paired:
        if len(n_sg) == 1:
            n = n_sg[0]
            n[n < 1e-20] = 1e-40
            sigma = sigma_xg[0]
            sigma[sigma < 1e-20] = 1e-40
            tau = tau_sg[0]
            tau[tau < 1e-20] = 1e-40

            # exchange
            e_x = twodexchange(n, sigma, tau, self.pars)
            e_g[:] += e_x * n

        # spin-polarized:
        else:
            n = n_sg
            n[n < 1e-20] = 1e-40
            sigma = sigma_xg
            sigma[sigma < 1e-20] = 1e-40
            tau = tau_sg
            tau[tau < 1e-20] = 1e-40

            # The spin polarized version is handle using the exact spin scaling
            # Ex[n1, n2] = (Ex[2*n1] + Ex[2*n2])/2
            na = 2.0 * n[0]
            nb = 2.0 * n[1]

            e2na_x = twodexchange(na, 4. * sigma[0], 2. * tau[0], self.pars)
            e2nb_x = twodexchange(nb, 4. * sigma[2], 2. * tau[1], self.pars)
            ea_x = e2na_x * na
            eb_x = e2nb_x * nb

            e_g[:] += (ea_x + eb_x) / 2.0


def twodexchange(n, sigma, tau, pars):
    # parameters for 2 Legendre polynomials
    parlen_i = pars[0]
    parlen_j = pars[2 + 2 * parlen_i]
    assert parlen_i == parlen_j
    pars_i = pars[1:2 + 2 * parlen_i]
    pars_j = pars[3 + 2 * parlen_i:]
    trans_i = pars_i[0]
    trans_j = pars_j[0]
    orders_i, coefs_i = np.split(pars_i[1:], 2)
    orders_j, coefs_j = np.split(pars_j[1:], 2)
    assert len(coefs_i) == len(orders_i)
    assert len(coefs_j) == len(orders_j)
    assert len(orders_i) == len(orders_j)

    # product Legendre expansion of Fx(s, alpha)
    e_x_ueg, rs = ueg_x(n)
    Fx = LegendreFx2(n, rs, sigma, tau,
                     trans_i, orders_i, coefs_i, trans_j, orders_j, coefs_j)
    return e_x_ueg * Fx


def LegendreFx2(n, rs, sigma, tau,
                trans_i, orders_i, coefs_i, trans_j, orders_j, coefs_j):
    # Legendre polynomial basis expansion in 2D

    # reduced density gradient in transformation t1(s)
    C2 = 0.26053088059892404
    s2 = sigma * (C2 * np.divide(rs, n))**2.
    x_i = transformation(s2, trans_i)
    assert(x_i.all() >= -1.0 and x_i.all() <= 1.0)

    # kinetic energy density parameter alpha in transformation t2(s)
    alpha = get_alpha(n, sigma, tau)
    x_j = transformation(alpha, trans_j)
    assert(x_j.all() >= -1.0 and x_j.all() <= 1.0)

    # product exchange enhancement factor
    Fx_i = legendre_polynomial(x_i, orders_i, coefs_i)
    Fx_j = legendre_polynomial(x_j, orders_j, coefs_j)
    Fx = Fx_i * Fx_j
    return Fx


def transformation(x, t):
    if t > 0:
        tmp = t + x
        x = 2.0 * np.divide(x, tmp) - 1.0
    elif int(t) == -1:
        tmp1 = (1.0 - x**2.0)**3.0
        tmp2 = (1.0 + x**3.0 + x**6.0)
        x = -1.0 * np.divide(tmp1, tmp2)
    else:
        raise KeyError('transformation %i unknown!' % t)
    return x


def get_alpha(n, sigma, tau):
    # tau LSDA
    aux = (3. / 10.) * (3.0 * np.pi * np.pi)**(2. / 3.)
    tau_lsda = aux * n**(5. / 3.)

    # von Weisaecker
    ind = (n != 0.).nonzero()
    gdms = np.maximum(sigma, 1e-40)  # |nabla rho|^2
    tau_w = np.zeros((np.shape(n)))
    tau_w[ind] = np.maximum(np.divide(gdms[ind], 8.0 * n[ind]), 1e-40)

    # z and alpha
    tau_ = np.maximum(tau_w, tau)
    z = np.divide(tau_w, tau_)
    alpha = np.divide(tau_ - tau_w, tau_lsda)
    assert(alpha.all() >= 0.0)
    return alpha


def ueg_x(n):
    C0I = 0.238732414637843
    C1 = -0.45816529328314287
    rs = (C0I / n)**(1 / 3.)
    ex = C1 / rs
    return ex, rs


def legendre_polynomial(x, orders, coefs, P=None):
    assert len(orders) == len(coefs)
    max_order = int(orders[-1])

    if P is None:
        P = np.zeros_like(x)
    else:
        assert np.shape(P) == np.shape(x)
    sh = np.shape(x)
    sh_ = np.append(sh, max_order + 2)
    L = np.empty(sh_)

    # initializing
    if len(sh) == 1:
        L[:, 0] = 1.0
        L[:, 1] = x
    else:
        L[:, :, :, 0] = 1.0
        L[:, :, :, 1] = x

    # recursively building polynomium terms
    if len(sh) == 1:
        for i in range(max_order):
            i += 2
            L[:, i] = (2.0 * x[:] * L[:, i - 1] - L[:, i - 2]
                - (x[:] * L[:, i - 1] - L[:, i - 2]) / i)
    else:
        for i in range(max_order):
            i += 2
            L[:, :, :, i] = (2.0 * x[:] * L[:, :, :, i - 1] - L[:, :, :, i - 2]
                - (x[:] * L[:, :, :, i - 1] - L[:, :, :, i - 2]) / i)

    # building polynomium P
    coefs_ = np.empty(max_order + 1)
    k = 0
    for i in range(len(coefs_)):
        if orders[k] == i:
            coefs_[i] = coefs[k]
            k += 1
        else:
            coefs_[i] = 0.0
    if len(sh) == 1:
        P += np.dot(L[:, :-1], coefs_)
    else:
        P += np.dot(L[:, :, :, :-1], coefs_)
    return P
