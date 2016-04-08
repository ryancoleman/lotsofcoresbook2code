import os
import sys
from time import time

import numpy as np
from scipy.special import sici
from scipy.special.orthogonal import p_roots

from ase.units import Hartree, Bohr
from ase.utils import prnt
from ase.utils.timing import timer

import gpaw.mpi as mpi
from gpaw.blacs import BlacsGrid, Redistributor
from gpaw.fd_operators import Gradient
from gpaw.io.tar import Writer, Reader
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.utilities.blas import gemmdot, axpy
from gpaw.wavefunctions.pw import PWDescriptor
from gpaw.xc.rpa import RPACorrelation


class FXCCorrelation(RPACorrelation):
    def __init__(self, calc, xc='RPA', filename=None,
                 skip_gamma=False, qsym=True, nlambda=8,
                 nfrequencies=16, frequency_max=800.0, frequency_scale=2.0,
                 frequencies=None, weights=None, density_cut=1.e-6,
                 world=mpi.world, nblocks=1,
                 unit_cells=None, tag=None,
                 txt=sys.stdout):

        RPACorrelation.__init__(self, calc, xc=xc, filename=filename,
                                skip_gamma=skip_gamma, qsym=qsym,
                                nfrequencies=nfrequencies, nlambda=nlambda,
                                frequency_max=frequency_max,
                                frequency_scale=frequency_scale,
                                frequencies=frequencies, weights=weights,
                                world=world, nblocks=nblocks,
                                txt=txt)

        self.l_l, self.weight_l = p_roots(nlambda)
        self.l_l = (self.l_l + 1.0) * 0.5
        self.weight_l *= 0.5
        self.xc = xc
        self.density_cut = density_cut
        if unit_cells is None:
            unit_cells = self.calc.wfs.kd.N_c
        self.unit_cells = unit_cells
        if tag is None:
            tag = self.calc.atoms.get_chemical_formula(mode='hill')
        self.tag = tag

    @timer('FXC')
    def calculate(self, ecut):
        if self.xc != 'RPA':
            if isinstance(ecut, (float, int)):
                self.ecut_max = ecut
            else:
                self.ecut_max = max(ecut)

            if not os.path.isfile('fhxc_%s_%s_%s_0.gpw'
                                  % (self.tag, self.xc, self.ecut_max)):
                kernel = Kernel(self.calc, self.xc, self.ibzq_qc,
                                self.fd, self.unit_cells, self.density_cut,
                                self.ecut_max, self.tag, self.timer)
                kernel.calculate_fhxc()
                del kernel
            else:
                prnt('%s kernel already calculated' % self.xc, file=self.fd)
                prnt(file=self.fd)

        if self.calc.wfs.nspins == 1:
            spin = False
        else:
            spin = True

        e = RPACorrelation.calculate(self, ecut, spin=spin)

        return e

    @timer('Chi0(q)')
    def calculate_q(self, chi0, pd,
                    chi0_swGG, chi0_swxvG, chi0_swvv,
                    Q_aGii, m1, m2, cut_G, A2_x):
        if chi0_swxvG is None:
            chi0_swxvG = range(2)  # Not used
            chi0_swvv = range(2)  # Not used
        chi0._calculate(pd, chi0_swGG[0], chi0_swxvG[0], chi0_swvv[0],
                        Q_aGii, m1, m2, [0])
        if len(chi0_swGG) == 2:
            chi0._calculate(pd, chi0_swGG[1], chi0_swxvG[1], chi0_swvv[1],
                            Q_aGii, m1, m2, [1])
        prnt('E_c(q) = ', end='', file=self.fd)

        chi0_swGG = chi0.redistribute(chi0_swGG, A2_x)
        
        if not pd.kd.gamma:
            e = self.calculate_energy(pd, chi0_swGG, cut_G)
            prnt('%.3f eV' % (e * Hartree), file=self.fd)
            self.fd.flush()
        else:
            e = 0.0
            for v in range(3):
                chi0_swGG[:, :, 0] = chi0_swxvG[:, :, 0, v]
                chi0_swGG[:, :, :, 0] = chi0_swxvG[:, :, 1, v]
                chi0_swGG[:, :, 0, 0] = chi0_swvv[:, :, v, v]
                ev = self.calculate_energy(pd, chi0_swGG, cut_G)
                e += ev
                prnt('%.3f' % (ev * Hartree), end='', file=self.fd)
                if v < 2:
                    prnt('/', end='', file=self.fd)
                else:
                    prnt('eV', file=self.fd)
                    self.fd.flush()
            e /= 3

        return e

    @timer('Energy')
    def calculate_energy(self, pd, chi0_swGG, cut_G):
        """Evaluate correlation energy from chi0 and the kernel fhxc"""

        ibzq2_q = [np.dot(self.ibzq_qc[i] - pd.kd.bzk_kc[0],
                          self.ibzq_qc[i] - pd.kd.bzk_kc[0])
                   for i in range(len(self.ibzq_qc))]
        qi = np.argsort(ibzq2_q)[0]

        G_G = pd.G2_qG[0]**0.5  # |G+q|

        if cut_G is not None:
            G_G = G_G[cut_G]

        nG = len(G_G)
        ns = len(chi0_swGG)

        if self.xc != 'RPA':
            r = Reader('fhxc_%s_%s_%s_%s.gpw' %
                       (self.tag, self.xc, self.ecut_max, qi))
            fv = r.get('fhxc_sGsG')
            if cut_G is not None:
                cut_sG = np.tile(cut_G, ns)
                cut_sG[len(cut_G):] += len(fv) / ns
                fv = fv.take(cut_sG, 0).take(cut_sG, 1)
            for s1 in range(ns):
                for s2 in range(ns):
                    m1 = s1 * nG
                    n1 = (s1 + 1) * nG
                    m2 = s2 * nG
                    n2 = (s2 + 1) * nG
                    fv[m1:n1, m2:n2] *= G_G * G_G[:, np.newaxis] / 4 / np.pi

                    if np.prod(self.unit_cells) > 1 and pd.kd.gamma:
                        m1 = s1 * nG
                        n1 = (s1 + 1) * nG
                        m2 = s2 * nG
                        n2 = (s2 + 1) * nG
                        fv[m1, m2:n2] = 0.0
                        fv[m1:n1, m2] = 0.0
                        fv[m1, m2] = 1.0
        else:
            fv = np.tile(np.eye(nG), (ns, ns))

        if pd.kd.gamma:
            G_G[0] = 1.0

        e_w = []
        for chi0_sGG in np.swapaxes(chi0_swGG, 0, 1):
            if cut_G is not None:
                chi0_sGG = chi0_sGG.take(cut_G, 1).take(cut_G, 2)
            chi0v = np.zeros((ns * nG, ns * nG), dtype=complex)
            for s in range(ns):
                m = s * nG
                n = (s + 1) * nG
                chi0v[m:n, m:n] = chi0_sGG[s] / G_G / G_G[:, np.newaxis]
            chi0v *= 4 * np.pi
            del chi0_sGG

            e = 0.0
            for l, weight in zip(self.l_l, self.weight_l):
                chiv = np.linalg.solve(np.eye(nG * ns) -
                                       l * np.dot(chi0v, fv), chi0v).real
                for s1 in range(ns):
                    for s2 in range(ns):
                        m1 = s1 * nG
                        n1 = (s1 + 1) * nG
                        m2 = s2 * nG
                        n2 = (s2 + 1) * nG
                        chiv_s1s2 = chiv[m1:n1, m2:n2]
                        e -= np.trace(chiv_s1s2) * weight
            e += np.trace(chi0v.real)
            e_w.append(e)
        E_w = np.zeros_like(self.omega_w)
        self.blockcomm.all_gather(np.array(e_w), E_w)
        energy = np.dot(E_w, self.weight_w) / (2 * np.pi)
        return energy


class Kernel:
    def __init__(self, calc, xc, ibzq_qc, fd, unit_cells,
                 density_cut, ecut, tag, timer):

        self.calc = calc
        self.gd = calc.density.gd
        self.xc = xc
        self.ibzq_qc = ibzq_qc
        self.fd = fd
        self.unit_cells = unit_cells
        self.density_cut = density_cut
        self.ecut = ecut
        self.tag = tag
        self.timer = timer
        
        self.A_x = -(3 / 4.) * (3 / np.pi)**(1 / 3.)

        self.n_g = calc.get_all_electron_density(gridrefinement=1)
        self.n_g *= Bohr**3

        if xc[-3:] == 'PBE':
            nf_g = calc.get_all_electron_density(gridrefinement=2)
            nf_g *= Bohr**3
            gdf = self.gd.refine()
            grad_v = [Gradient(gdf, v, n=1).apply for v in range(3)]
            gradnf_vg = gdf.empty(3)
            for v in range(3):
                grad_v[v](nf_g, gradnf_vg[v])
            self.gradn_vg = gradnf_vg[:, ::2, ::2, ::2]

        qd = KPointDescriptor(self.ibzq_qc)
        self.pd = PWDescriptor(ecut / Hartree, self.gd, complex, qd)

    @timer('FHXC')
    def calculate_fhxc(self):

        prnt('Calculating %s kernel at %d eV cutoff' %
             (self.xc, self.ecut), file=self.fd)
        if self.xc[0] == 'r':
            self.calculate_rkernel()
        else:
            assert self.xc[0] == 'A'
            self.calculate_local_kernel()

    def calculate_rkernel(self):

        gd = self.gd
        ng_c = gd.N_c
        cell_cv = gd.cell_cv
        icell_cv = 2 * np.pi * np.linalg.inv(cell_cv)
        vol = np.linalg.det(cell_cv)

        ns = self.calc.wfs.nspins
        n_g = self.n_g   # density on rough grid

        fx_g = ns * self.get_fxc_g(n_g)   # local exchange kernel
        qc_g = (-4 * np.pi * ns / fx_g)**0.5   # cutoff functional
        flocal_g = qc_g**3 * fx_g / (6 * np.pi**2)   # ren. x-kernel for r=r'
        Vlocal_g = 2 * qc_g / np.pi   # ren. Hartree kernel for r=r'

        ng = np.prod(ng_c)   # number of grid points
        r_vg = gd.get_grid_point_coordinates()
        rx_g = r_vg[0].flatten()
        ry_g = r_vg[1].flatten()
        rz_g = r_vg[2].flatten()

        prnt('    %d grid points and %d plane waves at the Gamma point' %
             (ng, self.pd.ngmax), file=self.fd)

        # Unit cells
        R_Rv = []
        weight_R = []
        nR_v = self.unit_cells
        nR = np.prod(nR_v)
        for i in range(-nR_v[0] + 1, nR_v[0]):
            for j in range(-nR_v[1] + 1, nR_v[1]):
                for h in range(-nR_v[2] + 1, nR_v[2]):
                    R_Rv.append(i * cell_cv[0] +
                                j * cell_cv[1] +
                                h * cell_cv[2])
                    weight_R.append((nR_v[0] - abs(i)) *
                                    (nR_v[1] - abs(j)) *
                                    (nR_v[2] - abs(h)) / float(nR))
        if nR > 1:
            # with more than one unit cell only the exchange kernel is
            # calculated on the grid. The bare Coulomb kernel is added
            # in PW basis and Vlocal_g only the exchange part
            dv = self.calc.density.gd.dv
            gc = (3 * dv / 4 / np.pi)**(1 / 3.)
            Vlocal_g -= 2 * np.pi * gc**2 / dv
            prnt('    Lattice point sampling: ' +
                 '(%s x %s x %s)^2 ' % (nR_v[0], nR_v[1], nR_v[2]) +
                 ' Reduced to %s lattice points' % len(R_Rv), file=self.fd)

        l_g_size = -(-ng // mpi.world.size)
        l_g_range = range(mpi.world.rank * l_g_size,
                          min((mpi.world.rank + 1) * l_g_size, ng))

        fhxc_qsGr = {}
        for iq in range(len(self.ibzq_qc)):
            fhxc_qsGr[iq] = np.zeros((ns, len(self.pd.G2_qG[iq]),
                                      len(l_g_range)), dtype=complex)

        inv_error = np.seterr()
        np.seterr(invalid='ignore')
        np.seterr(divide='ignore')

        t0 = time()
        # Loop over Lattice points
        for i, R_v in enumerate(R_Rv):
            # Loop over r'. f_rr and V_rr are functions of r (dim. as r_vg[0])
            if i == 1:
                prnt('      Finished 1 cell in %s seconds' % int(time() - t0) +
                     ' - estimated %s seconds left' %
                     int((len(R_Rv) - 1) * (time() - t0)),
                     file=self.fd)
                self.fd.flush()
            if len(R_Rv) > 5:
                if (i + 1) % (len(R_Rv) / 5 + 1) == 0:
                    prnt('      Finished %s cells in %s seconds'
                         % (i, int(time() - t0))
                         + ' - estimated %s seconds left'
                         % int((len(R_Rv) - i) * (time() - t0) / i),
                         file=self.fd)
                    self.fd.flush()
            for g in l_g_range:
                rx = rx_g[g] + R_v[0]
                ry = ry_g[g] + R_v[1]
                rz = rz_g[g] + R_v[2]

                # |r-r'-R_i|
                rr = ((r_vg[0] - rx)**2 +
                      (r_vg[1] - ry)**2 +
                      (r_vg[2] - rz)**2)**0.5

                n_av = (n_g + n_g.flatten()[g]) / 2.
                fx_g = ns * self.get_fxc_g(n_av, index=g)
                qc_g = (-4 * np.pi * ns / fx_g)**0.5
                x = qc_g * rr
                osc_x = np.sin(x) - x * np.cos(x)
                f_rr = fx_g * osc_x / (2 * np.pi**2 * rr**3)
                if nR > 1:   # include only exchange part of the kernel here
                    V_rr = (sici(x)[0] * 2 / np.pi - 1) / rr
                else:        # include the full kernel (also hartree part)
                    V_rr = (sici(x)[0] * 2 / np.pi) / rr

                # Terms with r = r'
                if (np.abs(R_v) < 0.001).all():
                    tmp_flat = f_rr.flatten()
                    tmp_flat[g] = flocal_g.flatten()[g]
                    f_rr = tmp_flat.reshape(ng_c)
                    tmp_flat = V_rr.flatten()
                    tmp_flat[g] = Vlocal_g.flatten()[g]
                    V_rr = tmp_flat.reshape(ng_c)
                    del tmp_flat

                f_rr[np.where(n_av < self.density_cut)] = 0.0
                V_rr[np.where(n_av < self.density_cut)] = 0.0

                f_rr *= weight_R[i]
                V_rr *= weight_R[i]

                # r-r'-R_i
                r_r = np.array([r_vg[0] - rx, r_vg[1] - ry, r_vg[2] - rz])

                # Fourier transform of r
                for iq, q in enumerate(self.ibzq_qc):
                    q_v = np.dot(q, icell_cv)
                    e_q = np.exp(-1j * gemmdot(q_v, r_r, beta=0.0))
                    f_q = self.pd.fft((f_rr + V_rr) * e_q, iq) * vol / ng
                    fhxc_qsGr[iq][0, :, g - l_g_range[0]] += f_q
                    if ns == 2:
                        f_q = self.pd.fft(V_rr * e_q, iq) * vol / ng
                        fhxc_qsGr[iq][1, :, g - l_g_range[0]] += f_q

        mpi.world.barrier()

        np.seterr(**inv_error)

        for iq, q in enumerate(self.ibzq_qc):
            npw = len(self.pd.G2_qG[iq])
            fhxc_sGsG = np.zeros((ns * npw, ns * npw), complex)
            l_pw_size = -(-npw // mpi.world.size)  # parallelize over PW below
            l_pw_range = range(mpi.world.rank * l_pw_size,
                               min((mpi.world.rank + 1) * l_pw_size, npw))

            if mpi.world.size > 1:
                # redistribute grid and plane waves in fhxc_qsGr[iq]
                bg1 = BlacsGrid(mpi.world, 1, mpi.world.size)
                bg2 = BlacsGrid(mpi.world, mpi.world.size, 1)
                bd1 = bg1.new_descriptor(npw, ng, npw, -(-ng / mpi.world.size))
                bd2 = bg2.new_descriptor(npw, ng, -(-npw / mpi.world.size), ng)

                fhxc_Glr = np.zeros((len(l_pw_range), ng), dtype=complex)
                if ns == 2:
                    Koff_Glr = np.zeros((len(l_pw_range), ng), dtype=complex)

                r = Redistributor(bg1.comm, bd1, bd2)
                r.redistribute(fhxc_qsGr[iq][0], fhxc_Glr, npw, ng)
                if ns == 2:
                    r.redistribute(fhxc_qsGr[iq][1], Koff_Glr, npw, ng)
            else:
                fhxc_Glr = fhxc_qsGr[iq][0]
                if ns == 2:
                    Koff_Glr = fhxc_qsGr[iq][1]

            # Fourier transform of r'
            for iG in range(len(l_pw_range)):
                f_g = fhxc_Glr[iG].reshape(ng_c)
                f_G = self.pd.fft(f_g.conj(), iq) * vol / ng
                fhxc_sGsG[l_pw_range[0] + iG, :npw] = f_G.conj()
                if ns == 2:
                    v_g = Koff_Glr[iG].reshape(ng_c)
                    v_G = self.pd.fft(v_g.conj(), iq) * vol / ng
                    fhxc_sGsG[npw + l_pw_range[0] + iG, :npw] = v_G.conj()

            if ns == 2:  # f_00 = f_11 and f_01 = f_10
                fhxc_sGsG[:npw, npw:] = fhxc_sGsG[npw:, :npw]
                fhxc_sGsG[npw:, npw:] = fhxc_sGsG[:npw, :npw]

            mpi.world.sum(fhxc_sGsG)
            fhxc_sGsG /= vol

            if mpi.rank == 0:
                w = Writer('fhxc_%s_%s_%s_%s.gpw' %
                           (self.tag, self.xc, self.ecut, iq))
                w.dimension('sG', ns * npw)
                w.add('fhxc_sGsG', ('sG', 'sG'), dtype=complex)
                if nR > 1:  # add Hartree kernel evaluated in PW basis
                    Gq2_G = self.pd.G2_qG[iq]
                    if (q == 0).all():
                        Gq2_G[0] = 1.
                    vq_G = 4 * np.pi / Gq2_G
                    fhxc_sGsG += np.tile(np.eye(npw) * vq_G, (ns, ns))
                w.fill(fhxc_sGsG)
                w.close()
            mpi.world.barrier()
        prnt(file=self.fd)

    def calculate_local_kernel(self):
        # Standard ALDA exchange kernel
        # Use with care. Results are very difficult to converge
        # Sensitive to density_cut
        ns = self.calc.wfs.nspins
        gd = self.gd
        pd = self.pd
        cell_cv = gd.cell_cv
        icell_cv = 2 * np.pi * np.linalg.inv(cell_cv)
        vol = np.linalg.det(cell_cv)

        fxc_sg = ns * self.get_fxc_g(ns * self.n_g)
        fxc_sg[np.where(self.n_g < self.density_cut)] = 0.0

        r_vg = gd.get_grid_point_coordinates()

        for iq in range(len(self.ibzq_qc)):
            Gvec_Gc = np.dot(pd.get_reciprocal_vectors(q=iq, add_q=False),
                             cell_cv / (2 * np.pi))
            npw = len(Gvec_Gc)
            l_pw_size = -(-npw // mpi.world.size)
            l_pw_range = range(mpi.world.rank * l_pw_size,
                               min((mpi.world.rank + 1) * l_pw_size, npw))
            fhxc_sGsG = np.zeros((ns * npw, ns * npw), dtype=complex)
            for s in range(ns):
                for iG in l_pw_range:
                    for jG in range(npw):
                        fxc = fxc_sg[s].copy()
                        dG_c = Gvec_Gc[iG] - Gvec_Gc[jG]
                        dG_v = np.dot(dG_c, icell_cv)
                        dGr_g = gemmdot(dG_v, r_vg, beta=0.0)
                        ft_fxc = gd.integrate(np.exp(-1j * dGr_g) * fxc)
                        fhxc_sGsG[s * npw + iG, s * npw + jG] = ft_fxc

            mpi.world.sum(fhxc_sGsG)
            fhxc_sGsG /= vol

            Gq2_G = self.pd.G2_qG[iq]
            if (self.ibzq_qc[iq] == 0).all():
                Gq2_G[0] = 1.
            vq_G = 4 * np.pi / Gq2_G
            fhxc_sGsG += np.tile(np.eye(npw) * vq_G, (ns, ns))

            if mpi.rank == 0:
                w = Writer('fhxc_%s_%s_%s_%s.gpw' %
                           (self.tag, self.xc, self.ecut, iq))
                w.dimension('sG', ns * npw)
                w.add('fhxc_sGsG', ('sG', 'sG'), dtype=complex)
                w.fill(fhxc_sGsG)
                w.close()
            mpi.world.barrier()
        prnt(file=self.fd)

    def get_fxc_g(self, n_g, index=None):
        if self.xc[-3:] == 'LDA':
            return self.get_lda_g(n_g)
        elif self.xc[-3:] == 'PBE':
            return self.get_pbe_g(n_g, index=index)
        else:
            raise '%s kernel not recognized' % self.xc

    def get_lda_g(self, n_g):
        return (4. / 9.) * self.A_x * n_g**(-2. / 3.)

    def get_pbe_g(self, n_g, index=None):
        if index is None:
            gradn_vg = self.gradn_vg
        else:
            gradn_vg = self.calc.density.gd.empty(3)
            for v in range(3):
                gradn_vg[v] = (self.gradn_vg[v] +
                               self.gradn_vg[v].flatten()[index]) / 2

        kf_g = (3. * np.pi**2 * n_g)**(1 / 3.)
        s2_g = np.zeros_like(n_g)
        for v in range(3):
            axpy(1.0, gradn_vg[v]**2, s2_g)
        s2_g /= 4 * kf_g**2 * n_g**2

        e_g = self.A_x * n_g**(4 / 3.)
        v_g = (4 / 3.) * e_g / n_g
        f_g = (1 / 3.) * v_g / n_g

        kappa = 0.804
        mu = 0.2195149727645171

        denom_g = (1 + mu * s2_g / kappa)
        F_g = 1. + kappa - kappa / denom_g
        Fn_g = -mu / denom_g**2 * 8 * s2_g / (3 * n_g)
        Fnn_g = -11 * Fn_g / (3 * n_g) - 2 * Fn_g**2 / kappa

        fxc_g = f_g * F_g
        fxc_g += 2 * v_g * Fn_g
        fxc_g += e_g * Fnn_g

        # Contributions from varying the gradient
        #Fgrad_vg = np.zeros_like(gradn_vg)
        #Fngrad_vg = np.zeros_like(gradn_vg)
        #for v in range(3):
        #    axpy(1.0, mu / den_g**2 * gradn_vg[v] / (2 * kf_g**2 * n_g**2),
        #         Fgrad_vg[v])
        #    axpy(-8.0, Fgrad_vg[v] / (3 * n_g), Fngrad_vg[v])
        #    axpy(-2.0, Fgrad_vg[v] * Fn_g / kappa, Fngrad_vg[v])

        #tmp = np.zeros_like(fxc_g)
        #tmp1 = np.zeros_like(fxc_g)

        #for v in range(3):
            #self.grad_v[v](Fgrad_vg[v], tmp)
            #axpy(-2.0, tmp * v_g, fxc_g)
            #for u in range(3):
                #self.grad_v[u](Fgrad_vg[u] * tmp, tmp1)
                #axpy(-4.0/kappa, tmp1 * e_g, fxc_g)
            #self.grad_v[v](Fngrad_vg[v], tmp)
            #axpy(-2.0, tmp * e_g, fxc_g)
        #self.laplace(mu / den_g**2 / (2 * kf_g**2 * n_g**2), tmp)
        #axpy(1.0, tmp * e_g, fxc_g)

        return fxc_g

"""
    def get_fxc_libxc_g(self, n_g):
        ### NOT USED AT THE MOMENT
        gd = self.calc.density.gd.refine()

        xc = XC('GGA_X_' + self.xc[2:])
        #xc = XC('LDA_X')
        #sigma = np.zeros_like(n_g).flat[:]
        xc.set_grid_descriptor(gd)
        sigma_xg, gradn_svg = xc.calculate_sigma(np.array([n_g]))

        dedsigma_xg = np.zeros_like(sigma_xg)
        e_g = np.zeros_like(n_g)
        v_sg = np.array([np.zeros_like(n_g)])

        xc.calculate_gga(e_g, np.array([n_g]), v_sg, sigma_xg, dedsigma_xg)

        sigma = sigma_xg[0].flat[:]
        gradn_vg = gradn_svg[0]
        dedsigma_g = dedsigma_xg[0]

        libxc = LibXC('GGA_X_' + self.xc[2:])
        #libxc = LibXC('LDA_X')
        libxc.initialize(1)
        libxc_fxc = libxc.xc.calculate_fxc_spinpaired

        fxc_g = np.zeros_like(n_g).flat[:]
        d2edndsigma_g = np.zeros_like(n_g).flat[:]
        d2ed2sigma_g = np.zeros_like(n_g).flat[:]

        libxc_fxc(n_g.flat[:], fxc_g, sigma, d2edndsigma_g, d2ed2sigma_g)
        fxc_g = fxc_g.reshape(np.shape(n_g))
        d2edndsigma_g = d2edndsigma_g.reshape(np.shape(n_g))
        d2ed2sigma_g = d2ed2sigma_g.reshape(np.shape(n_g))

        tmp = np.zeros_like(fxc_g)
        tmp1 = np.zeros_like(fxc_g)

        #for v in range(3):
            #self.grad_v[v](d2edndsigma_g * gradn_vg[v], tmp)
            #axpy(-4.0, tmp, fxc_g)

        #for u in range(3):
            #for v in range(3):
                #self.grad_v[v](d2ed2sigma_g * gradn_vg[u] * gradn_vg[v], tmp)
                #self.grad_v[u](tmp, tmp1)
                #axpy(4.0, tmp1, fxc_g)

        #self.laplace(dedsigma_g, tmp)
        #axpy(2.0, tmp, fxc_g)

        return fxc_g[::2, ::2, ::2]

    def get_numerical_fxc_sg(self, n_sg):
        ### NOT USED AT THE MOMENT
        gd = self.calc.density.gd.refine()
        delta = 1.e-4

        if self.xc[2:] == 'LDA':
            xc = XC('LDA_X')
            v1xc_sg = np.zeros_like(n_sg)
            v2xc_sg = np.zeros_like(n_sg)
            xc.calculate(gd, (1 + delta) * n_sg, v1xc_sg)
            xc.calculate(gd, (1 - delta) * n_sg, v2xc_sg)
            fxc_sg = (v1xc_sg - v2xc_sg) / (2 * delta * n_sg)
        else:
            fxc_sg = np.zeros_like(n_sg)
            xc = XC('GGA_X_' + self.xc[2:])
            vxc_sg = np.zeros_like(n_sg)
            xc.calculate(gd, n_sg, vxc_sg)
            for s in range(len(n_sg)):
                for x in range(len(n_sg[0])):
                    for y in range(len(n_sg[0, 0])):
                        for z in range(len(n_sg[0, 0, 0])):
                            v1xc_sg = np.zeros_like(n_sg)
                            n1_sg = n_sg.copy()
                            n1_sg[s, x, y, z] *= (1 + delta)
                            xc.calculate(gd, n1_sg, v1xc_sg)
                            num = v1xc_sg[s, x, y, z] - vxc_sg[s, x, y, z]
                            den = delta * n_sg[s, x, y, z]
                            fxc_sg[s, x, y, z] = num / den

        return fxc_sg[:, ::2, ::2, ::2]
"""
