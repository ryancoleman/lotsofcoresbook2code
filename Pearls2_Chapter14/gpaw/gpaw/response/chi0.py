from __future__ import print_function, division

import sys
from time import ctime

import numpy as np
from ase.units import Hartree
from ase.utils.timing import timer

import gpaw.mpi as mpi
from gpaw import extra_parameters
from gpaw.blacs import (BlacsGrid, BlacsDescriptor, Redistributor,
                        DryRunBlacsGrid)
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.occupations import FermiDirac
from gpaw.response.pair import PairDensity
from gpaw.utilities.memory import maxrss
from gpaw.utilities.progressbar import ProgressBar
from gpaw.utilities.blas import gemm, rk, czher, mmm
from gpaw.wavefunctions.pw import PWDescriptor


def frequency_grid(domega0, omega2, omegamax):
    beta = (2**0.5 - 1) * domega0 / omega2
    wmax = int(omegamax / (domega0 + beta * omegamax)) + 2
    w = np.arange(wmax)
    omega_w = w * domega0 / (1 - beta * w)
    return omega_w


class Chi0(PairDensity):
    def __init__(self, calc,
                 frequencies=None, domega0=0.1, omega2=10.0, omegamax=None,
                 ecut=50, hilbert=True, nbands=None,
                 timeordered=False, eta=0.2, ftol=1e-6, threshold=1,
                 real_space_derivatives=False, intraband=True,
                 world=mpi.world, txt=sys.stdout, timer=None,
                 nblocks=1, no_optical_limit=False,
                 keep_occupied_states=False, gate_voltage=None):

        PairDensity.__init__(self, calc, ecut, ftol, threshold,
                             real_space_derivatives, world, txt, timer,
                             nblocks=nblocks,
                             gate_voltage=gate_voltage)

        self.eta = eta / Hartree
        self.domega0 = domega0 / Hartree
        self.omega2 = omega2 / Hartree
        self.omegamax = None if omegamax is None else omegamax / Hartree
        self.nbands = nbands or self.calc.wfs.bd.nbands
        self.keep_occupied_states = keep_occupied_states
        self.intraband = intraband
        self.no_optical_limit = no_optical_limit
        
        omax = self.find_maximum_frequency()

        if frequencies is None:
            if self.omegamax is None:
                self.omegamax = omax
            print('Using nonlinear frequency grid from 0 to %.3f eV' %
                  (self.omegamax * Hartree), file=self.fd)
            self.omega_w = frequency_grid(self.domega0, self.omega2,
                                          self.omegamax)
        else:
            self.omega_w = np.asarray(frequencies) / Hartree
            assert not hilbert
        self.hilbert = hilbert
        self.timeordered = bool(timeordered)

        if self.eta == 0.0:
            assert not hilbert
            assert not timeordered
            assert not self.omega_w.real.any()

        # Occupied states:
        self.mysKn1n2 = None  # my (s, K, n1, n2) indices
        self.distribute_k_points_and_bands(0, self.nocc2)
        self.mykpts = None

        wfs = self.calc.wfs
        self.prefactor = 2 / self.vol / wfs.kd.nbzkpts / wfs.nspins

        self.chi0_vv = None  # strength of intraband peak

    def find_maximum_frequency(self):
        self.epsmin = 10000.0
        self.epsmax = -10000.0
        for kpt in self.calc.wfs.kpt_u:
            self.epsmin = min(self.epsmin, kpt.eps_n[0])
            self.epsmax = max(self.epsmax, kpt.eps_n[self.nbands - 1])

        print('Minimum eigenvalue: %10.3f eV' % (self.epsmin * Hartree),
              file=self.fd)
        print('Maximum eigenvalue: %10.3f eV' % (self.epsmax * Hartree),
              file=self.fd)

        return self.epsmax - self.epsmin

    def calculate(self, q_c, spin='all', A_x=None):
        wfs = self.calc.wfs

        if spin == 'all':
            spins = range(wfs.nspins)
        else:
            assert spin in range(wfs.nspins)
            spins = [spin]

        q_c = np.asarray(q_c, dtype=float)
        qd = KPointDescriptor([q_c])
        pd = PWDescriptor(self.ecut, wfs.gd, complex, qd)

        self.print_chi(pd)

        if extra_parameters.get('df_dry_run'):
            print('    Dry run exit', file=self.fd)
            raise SystemExit

        nG = pd.ngmax
        nw = len(self.omega_w)
        mynG = (nG + self.blockcomm.size - 1) // self.blockcomm.size
        self.Ga = self.blockcomm.rank * mynG
        self.Gb = min(self.Ga + mynG, nG)
        assert mynG * (self.blockcomm.size - 1) < nG
        
        if A_x is not None:
            nx = nw * (self.Gb - self.Ga) * nG
            chi0_wGG = A_x[:nx].reshape((nw, self.Gb - self.Ga, nG))
            chi0_wGG[:] = 0.0
        else:
            chi0_wGG = np.zeros((nw, self.Gb - self.Ga, nG), complex)

        if np.allclose(q_c, 0.0):
            chi0_wxvG = np.zeros((len(self.omega_w), 2, 3, nG), complex)
            chi0_wvv = np.zeros((len(self.omega_w), 3, 3), complex)
            self.chi0_vv = np.zeros((3, 3), complex)
        else:
            chi0_wxvG = None
            chi0_wvv = None

        print('Initializing PAW Corrections', file=self.fd)
        self.Q_aGii = self.initialize_paw_corrections(pd)

        # Do all empty bands:
        m1 = self.nocc1
        m2 = self.nbands
        
        self._calculate(pd, chi0_wGG, chi0_wxvG, chi0_wvv, self.Q_aGii,
                        m1, m2, spins)
        
        return pd, chi0_wGG, chi0_wxvG, chi0_wvv

    @timer('Calculate CHI_0')
    def _calculate(self, pd, chi0_wGG, chi0_wxvG, chi0_wvv, Q_aGii,
                   m1, m2, spins):
        wfs = self.calc.wfs

        if self.keep_occupied_states:
            self.mykpts = [self.get_k_point(s, K, n1, n2)
                           for s, K, n1, n2 in self.mysKn1n2]

        if self.eta == 0.0:
            update = self.update_hermitian
        elif self.hilbert:
            update = self.update_hilbert
        else:
            update = self.update

        q_c = pd.kd.bzk_kc[0]
        optical_limit = not self.no_optical_limit and np.allclose(q_c, 0.0)

        pb = ProgressBar(self.fd)

        self.timer.start('Loop')
        # kpt1 occupied and kpt2 empty:
        for kn, (s, K, n1, n2) in enumerate(self.mysKn1n2):
            pb.update(kn / len(self.mysKn1n2))
            if self.keep_occupied_states:
                kpt1 = self.mykpts[kn]
            else:
                kpt1 = self.get_k_point(s, K, n1, n2)

            if kpt1.s not in spins:
                continue

            with self.timer('k+q'):
                K2 = wfs.kd.find_k_plus_q(q_c, [kpt1.K])[0]
            with self.timer('get k2'):
                kpt2 = self.get_k_point(kpt1.s, K2, m1, m2, block=True)
            with self.timer('fft-indices'):
                Q_G = self.get_fft_indices(kpt1.K, kpt2.K, q_c, pd,
                                           kpt1.shift_c - kpt2.shift_c)

            for n in range(kpt1.n2 - kpt1.n1):
                eps1 = kpt1.eps_n[n]
                
                # Only update if there exists deps <= omegamax
                if self.omegamax is not None:
                    m = [m for m, d in enumerate(eps1 - kpt2.eps_n)
                         if abs(d) <= self.omegamax]
                else:
                    m = range(0, kpt2.n2 - kpt2.n1)
                
                if not len(m):
                    continue

                deps_m = (eps1 - kpt2.eps_n)[m]
                f1 = kpt1.f_n[n]
                with self.timer('conj'):
                    ut1cc_R = kpt1.ut_nR[n].conj()
                with self.timer('paw'):
                    C1_aGi = [np.dot(Q_Gii, P1_ni[n].conj())
                              for Q_Gii, P1_ni in zip(Q_aGii, kpt1.P_ani)]
                n_mG = self.calculate_pair_densities(ut1cc_R, C1_aGi, kpt2,
                                                     pd, Q_G)[m]
                df_m = (f1 - kpt2.f_n)[m]

                # This is not quite right for degenerate partially occupied
                # bands, but good enough for now:
                df_m[df_m <= 1e-20] = 0.0

                if optical_limit:
                    self.update_optical_limit(
                        n, m, kpt1, kpt2, deps_m, df_m, n_mG,
                        chi0_wxvG, chi0_wvv)

                update(n_mG, deps_m, df_m, chi0_wGG)

            if optical_limit and self.intraband:
                # Avoid that more ranks are summing up
                # the intraband contributions
                if kpt1.n1 == 0 and self.blockcomm.rank == 0:
                    assert self.nocc2 <= kpt2.nb, \
                        print('Error: Too few unoccupied bands')
                    self.update_intraband(kpt2)

        self.timer.stop('Loop')

        pb.finish()

        with self.timer('Sum CHI_0'):
            for chi0_GG in chi0_wGG:
                self.kncomm.sum(chi0_GG)

            if optical_limit:
                self.kncomm.sum(chi0_wxvG)
                self.kncomm.sum(chi0_wvv)
                if self.intraband:
                    self.kncomm.sum(self.chi0_vv)

        print('Memory used: {0:.3f} MB / CPU'.format(maxrss() / 1024**2),
              file=self.fd)

        if (self.eta == 0.0 or self.hilbert) and self.blockcomm.size == 1:
            # Fill in upper/lower triangle also:
            nG = pd.ngmax
            il = np.tril_indices(nG, -1)
            iu = il[::-1]
            if self.hilbert:
                for chi0_GG in chi0_wGG:
                    chi0_GG[il] = chi0_GG[iu].conj()
            else:
                for chi0_GG in chi0_wGG:
                    chi0_GG[iu] = chi0_GG[il].conj()

        if self.hilbert:
            with self.timer('Hilbert transform'):
                ht = HilbertTransform(self.omega_w, self.eta,
                                      self.timeordered)
                ht(chi0_wGG)
                if optical_limit:
                    ht(chi0_wvv)
                    ht(chi0_wxvG)
            print('Hilbert transform done', file=self.fd)

        if optical_limit and self.intraband:  # Add intraband contribution
            omega_w = self.omega_w.copy()
            if omega_w[0] == 0.0:
                omega_w[0] = 1e-14

            chi0_vv = self.chi0_vv
            self.world.broadcast(chi0_vv, 0)

            chi0_wvv += (chi0_vv[np.newaxis] /
                         (omega_w[:, np.newaxis, np.newaxis] *
                          (omega_w[:, np.newaxis, np.newaxis] +
                           1j * self.eta)))

        return pd, chi0_wGG, chi0_wxvG, chi0_wvv

    @timer('CHI_0 update')
    def update(self, n_mG, deps_m, df_m, chi0_wGG):
        if self.timeordered:
            deps1_m = deps_m + 1j * self.eta * np.sign(deps_m)
            deps2_m = deps1_m
        else:
            deps1_m = deps_m + 1j * self.eta
            deps2_m = deps_m - 1j * self.eta

        for omega, chi0_GG in zip(self.omega_w, chi0_wGG):
            x_m = df_m * (1 / (omega + deps1_m) - 1 / (omega - deps2_m))
            nx_mG = n_mG * x_m[:, np.newaxis]
            gemm(self.prefactor, n_mG.conj(), np.ascontiguousarray(nx_mG.T),
                 1.0, chi0_GG)

    @timer('CHI_0 hermetian update')
    def update_hermitian(self, n_mG, deps_m, df_m, chi0_wGG):
        for w, omega in enumerate(self.omega_w):
            if self.blockcomm.size == 1:
                x_m = (-2 * df_m * deps_m / (omega.imag**2 + deps_m**2))**0.5
                nx_mG = n_mG.conj() * x_m[:, np.newaxis]
                rk(-self.prefactor, nx_mG, 1.0, chi0_wGG[w], 'n')
            else:
                x_m = 2 * df_m * deps_m / (omega.imag**2 + deps_m**2)
                mynx_mG = n_mG[:, self.Ga:self.Gb] * x_m[:, np.newaxis]
                mmm(self.prefactor, mynx_mG, 'c', n_mG, 'n', 1.0, chi0_wGG[w])

    @timer('CHI_0 spectral function update')
    def update_hilbert(self, n_mG, deps_m, df_m, chi0_wGG):
        self.timer.start('prep')
        beta = (2**0.5 - 1) * self.domega0 / self.omega2
        o_m = abs(deps_m)
        w_m = (o_m / (self.domega0 + beta * o_m)).astype(int)
        o1_m = self.omega_w[w_m]
        o2_m = self.omega_w[w_m + 1]
        p_m = self.prefactor * abs(df_m) / (o2_m - o1_m)**2  # XXX abs()?
        p1_m = p_m * (o2_m - o_m)
        p2_m = p_m * (o_m - o1_m)
        self.timer.stop('prep')

        if self.blockcomm.size > 1:
            for p1, p2, n_G, w in zip(p1_m, p2_m, n_mG, w_m):
                myn_G = n_G[self.Ga:self.Gb].reshape((-1, 1))
                gemm(p1, n_G.reshape((-1, 1)), myn_G, 1.0, chi0_wGG[w], 'c')
                gemm(p2, n_G.reshape((-1, 1)), myn_G, 1.0, chi0_wGG[w + 1],
                     'c')
            return

        for p1, p2, n_G, w in zip(p1_m, p2_m, n_mG, w_m):
            czher(p1, n_G.conj(), chi0_wGG[w])
            czher(p2, n_G.conj(), chi0_wGG[w + 1])

    @timer('CHI_0 optical limit update')
    def update_optical_limit(self, n, m, kpt1, kpt2, deps_m, df_m, n_mG,
                             chi0_wxvG, chi0_wvv):
        n0_mv = PairDensity.update_optical_limit(self, n, m, kpt1, kpt2,
                                                 deps_m, df_m, n_mG)

        if self.hilbert:
            self.update_optical_limit_hilbert(n0_mv, deps_m, df_m, n_mG,
                                              chi0_wxvG, chi0_wvv)
            return

        if self.timeordered:
            # avoid getting a zero from np.sign():
            deps1_m = deps_m + 1j * self.eta * np.sign(deps_m + 1e-20)
            deps2_m = deps1_m
        else:
            deps1_m = deps_m + 1j * self.eta
            deps2_m = deps_m - 1j * self.eta

        for w, omega in enumerate(self.omega_w):
            x_m = self.prefactor * df_m * (1 / (omega + deps1_m) -
                                           1 / (omega - deps2_m))

            chi0_wvv[w] += np.dot(x_m * n0_mv.T, n0_mv.conj())
            chi0_wxvG[w, 0, :, 1:] += np.dot(x_m * n0_mv.T, n_mG[:, 1:].conj())
            chi0_wxvG[w, 1, :, 1:] += np.dot(x_m * n0_mv.T.conj(), n_mG[:, 1:])

    @timer('CHI_0 optical limit hilbert-update')
    def update_optical_limit_hilbert(self, n0_mv, deps_m, df_m, n_mG,
                                     chi0_wxvG, chi0_wvv):
        beta = (2**0.5 - 1) * self.domega0 / self.omega2
        for deps, df, n0_v, n_G in zip(deps_m, df_m, n0_mv, n_mG):
            o = abs(deps)
            w = int(o / (self.domega0 + beta * o))
            if w + 2 > len(self.omega_w):
                break
            o1, o2 = self.omega_w[w:w + 2]
            assert o1 <= o <= o2, (o1, o, o2)

            p = self.prefactor * abs(df) / (o2 - o1)**2  # XXX abs()?
            p1 = p * (o2 - o)
            p2 = p * (o - o1)
            x_vv = np.outer(n0_v, n0_v.conj())
            chi0_wvv[w] += p1 * x_vv
            chi0_wvv[w + 1] += p2 * x_vv
            x_vG = np.outer(n0_v, n_G[1:].conj())
            chi0_wxvG[w, 0, :, 1:] += p1 * x_vG
            chi0_wxvG[w + 1, 0, :, 1:] += p2 * x_vG
            chi0_wxvG[w, 1, :, 1:] += p1 * x_vG.conj()
            chi0_wxvG[w + 1, 1, :, 1:] += p2 * x_vG.conj()

    @timer('CHI_0 intraband update')
    def update_intraband(self, kpt):
        """Check whether there are any partly occupied bands."""

        width = self.calc.occupations.width
        if width == 0.0:
            return

        assert isinstance(self.calc.occupations, FermiDirac)
        dfde_m = - 1. / width * (kpt.f_n - kpt.f_n**2.0)
        partocc_m = np.abs(dfde_m) > 1e-5
        if not partocc_m.any():
            return

        # Break bands into degenerate chunks
        deginds_cm = []  # indexing c as chunk number
        for m in range(kpt.nb - kpt.na):
            inds_m = np.nonzero(np.abs(kpt.eps_n[m] - kpt.eps_n) < 1e-5)[0]
            if m == np.min(inds_m) and partocc_m[m]:
                deginds_cm.append((inds_m))

        # Sum over the chunks of degenerate bands
        for inds_m in deginds_cm:
            deg = len(inds_m)
            vel_mmv = -1j * PairDensity.update_intraband(self, inds_m, kpt)
            vel_mv = np.zeros((deg, 3), dtype=complex)

            for iv in range(3):
                w, v = np.linalg.eig(vel_mmv[..., iv])
                vel_mv[:, iv] = w

            for m in range(deg):
                velm_v = vel_mv[m]
                x_vv = (-self.prefactor * dfde_m[inds_m[m]] *
                        np.outer(velm_v.conj(), velm_v))

                self.chi0_vv += x_vv

    @timer('redist')
    def redistribute(self, in_wGG, out_x=None):
        """Redistribute array.
        
        Switch between two kinds of parallel distributions:
            
        1) parallel over G-vectors (second dimension of in_wGG)
        2) parallel over frequency (first dimension of in_wGG)

        Returns new array using the memory in the 1-d array out_x.
        """
        
        comm = self.blockcomm
        
        if comm.size == 1:
            return in_wGG
            
        nw = len(self.omega_w)
        nG = in_wGG.shape[2]
        mynw = (nw + comm.size - 1) // comm.size
        mynG = (nG + comm.size - 1) // comm.size
        
        bg1 = BlacsGrid(comm, comm.size, 1)
        bg2 = BlacsGrid(comm, 1, comm.size)
        md1 = BlacsDescriptor(bg1, nw, nG**2, mynw, nG**2)
        md2 = BlacsDescriptor(bg2, nw, nG**2, nw, mynG * nG)
        
        if len(in_wGG) == nw:
            mdin = md2
            mdout = md1
        else:
            mdin = md1
            mdout = md2
            
        r = Redistributor(comm, mdin, mdout)
        
        outshape = (mdout.shape[0], mdout.shape[1] // nG, nG)
        if out_x is None:
            out_wGG = np.empty(outshape, complex)
        else:
            out_wGG = out_x[:np.product(outshape)].reshape(outshape)

        r.redistribute(in_wGG.reshape(mdin.shape),
                       out_wGG.reshape(mdout.shape))
        
        return out_wGG

    @timer('dist freq')
    def distribute_frequencies(self, chi0_wGG):
        """Distribute frequencies to all cores."""
        
        world = self.world
        comm = self.blockcomm
        
        if world.size == 1:
            return chi0_wGG
            
        nw = len(self.omega_w)
        nG = chi0_wGG.shape[2]
        mynw = (nw + world.size - 1) // world.size
        mynG = (nG + comm.size - 1) // comm.size
  
        wa = min(world.rank * mynw, nw)
        wb = min(wa + mynw, nw)

        if self.blockcomm.size == 1:
            return chi0_wGG[wa:wb].copy()

        if self.kncomm.rank == 0:
            bg1 = BlacsGrid(comm, 1, comm.size)
            in_wGG = chi0_wGG.reshape((nw, -1))
        else:
            bg1 = DryRunBlacsGrid(mpi.serial_comm, 1, 1)
            in_wGG = np.zeros((0, 0), complex)
        md1 = BlacsDescriptor(bg1, nw, nG**2, nw, mynG * nG)
        
        bg2 = BlacsGrid(world, world.size, 1)
        md2 = BlacsDescriptor(bg2, nw, nG**2, mynw, nG**2)
        
        r = Redistributor(world, md1, md2)
        shape = (wb - wa, nG, nG)
        out_wGG = np.empty(shape, complex)
        r.redistribute(in_wGG, out_wGG.reshape((wb - wa, nG**2)))
        
        return out_wGG

    def print_chi(self, pd):
        calc = self.calc
        gd = calc.wfs.gd

        if extra_parameters.get('df_dry_run'):
            from gpaw.mpi import DryRunCommunicator
            size = extra_parameters['df_dry_run']
            world = DryRunCommunicator(size)
        else:
            world = self.world

        print('%s' % ctime(), file=self.fd)
        print('Called response.chi0.calculate with', file=self.fd)

        q_c = pd.kd.bzk_kc[0]
        print('    q_c: [%f, %f, %f]' % (q_c[0], q_c[1], q_c[2]), file=self.fd)

        nw = len(self.omega_w)
        print('    Number of frequency points: %d' % nw, file=self.fd)

        ecut = self.ecut * Hartree
        print('    Planewave cutoff: %f' % ecut, file=self.fd)

        ns = calc.wfs.nspins
        print('    Number of spins: %d' % ns, file=self.fd)

        nbands = self.nbands
        print('    Number of bands: %d' % nbands, file=self.fd)

        nk = calc.wfs.kd.nbzkpts
        print('    Number of kpoints: %d' % nk, file=self.fd)

        nik = calc.wfs.kd.nibzkpts
        print('    Number of irredicible kpoints: %d' % nik, file=self.fd)
        
        ngmax = pd.ngmax
        print('    Number of planewaves: %d' % ngmax, file=self.fd)

        eta = self.eta * Hartree
        print('    Broadening (eta): %f' % eta, file=self.fd)
        
        wsize = world.size
        print('    world.size: %d' % wsize, file=self.fd)

        knsize = self.kncomm.size
        print('    kncomm.size: %d' % knsize, file=self.fd)

        bsize = self.blockcomm.size
        print('    blockcomm.size: %d' % bsize, file=self.fd)
        
        nocc = self.nocc1
        print('    Number of completely occupied states: %d'
              % nocc, file=self.fd)
        
        npocc = self.nocc2
        print('    Number of partially occupied states: %d'
              % npocc, file=self.fd)

        keep = self.keep_occupied_states
        print('    Keep occupied states: %s' % keep, file=self.fd)

        print('', file=self.fd)
        print('    Memory estimate of potentially large arrays:', file=self.fd)

        chisize = nw * pd.ngmax**2 * 16. / 1024**2
        print('        chi0_wGG: %f M / cpu' % chisize, file=self.fd)

        ngridpoints = gd.N_c[0] * gd.N_c[1] * gd.N_c[2]

        if self.keep_occupied_states:
            nstat = (ns * nk * npocc + world.size - 1) // world.size
        else:
            nstat = (ns * npocc + world.size - 1) // world.size

        occsize = nstat * ngridpoints * 16. / 1024**2
        print('        Occupied states: %f M / cpu' % occsize,
              file=self.fd)

        print('        Memory usage before allocation: %f M / cpu'
              % (maxrss() / 1024**2), file=self.fd)

        print('', file=self.fd)


class HilbertTransform:
    def __init__(self, omega_w, eta, timeordered=False, gw=False,
                 blocksize=500):
        """Analytic Hilbert transformation using linear interpolation.

        Hilbert transform::

           oo
          /           1                1
          |dw' (-------------- - --------------) S(w').
          /     w - w' + i eta   w + w' + i eta
          0

        With timeordered=True, you get::

           oo
          /           1                1
          |dw' (-------------- - --------------) S(w').
          /     w - w' - i eta   w + w' + i eta
          0

        With gw=True, you get::

           oo
          /           1                1
          |dw' (-------------- + --------------) S(w').
          /     w - w' + i eta   w + w' + i eta
          0

        """

        self.blocksize = blocksize

        if timeordered:
            self.H_ww = self.H(omega_w, -eta) + self.H(omega_w, -eta, -1)
        elif gw:
            self.H_ww = self.H(omega_w, eta) - self.H(omega_w, -eta, -1)
        else:
            self.H_ww = self.H(omega_w, eta) + self.H(omega_w, -eta, -1)

    def H(self, o_w, eta, sign=1):
        """Calculate transformation matrix.

        With s=sign (+1 or -1)::

                        oo
                       /       dw'
          X (w, eta) = | ---------------- S(w').
           s           / s w - w' + i eta
                       0

        Returns H_ij so that X_i = np.dot(H_ij, S_j), where::

            X_i = X (omega_w[i]) and S_j = S(omega_w[j])
                   s
        """

        nw = len(o_w)
        H_ij = np.zeros((nw, nw), complex)
        do_j = o_w[1:] - o_w[:-1]
        for i, o in enumerate(o_w):
            d_j = o_w - o * sign
            y_j = 1j * np.arctan(d_j / eta) + 0.5 * np.log(d_j**2 + eta**2)
            y_j = (y_j[1:] - y_j[:-1]) / do_j
            H_ij[i, :-1] = 1 - (d_j[1:] - 1j * eta) * y_j
            H_ij[i, 1:] -= 1 - (d_j[:-1] - 1j * eta) * y_j
        return H_ij

    def __call__(self, S_wx):
        """Inplace transform"""
        B_wx = S_wx.reshape((len(S_wx), -1))
        nw, nx = B_wx.shape
        tmp_wx = np.zeros((nw, min(nx, self.blocksize)), complex)
        for x in range(0, nx, self.blocksize):
            b_wx = B_wx[:, x:x + self.blocksize]
            c_wx = tmp_wx[:, :b_wx.shape[1]]
            gemm(1.0, b_wx, self.H_ww, 0.0, c_wx)
            b_wx[:] = c_wx


if __name__ == '__main__':
    do = 0.025
    eta = 0.1
    omega_w = frequency_grid(do, 10.0, 3)
    print(len(omega_w))
    X_w = omega_w * 0j
    Xt_w = omega_w * 0j
    Xh_w = omega_w * 0j
    for o in -np.linspace(2.5, 2.9, 10):
        X_w += (1 / (omega_w + o + 1j * eta) -
                1 / (omega_w - o + 1j * eta)) / o**2
        Xt_w += (1 / (omega_w + o - 1j * eta) -
                 1 / (omega_w - o + 1j * eta)) / o**2
        w = int(-o / do / (1 + 3 * -o / 10))
        o1, o2 = omega_w[w:w + 2]
        assert o1 - 1e-12 <= -o <= o2 + 1e-12, (o1, -o, o2)
        p = 1 / (o2 - o1)**2 / o**2
        Xh_w[w] += p * (o2 - -o)
        Xh_w[w + 1] += p * (-o - o1)

    ht = HilbertTransform(omega_w, eta, 1)
    ht(Xh_w)

    import matplotlib.pyplot as plt
    plt.plot(omega_w, X_w.imag, label='ImX')
    plt.plot(omega_w, X_w.real, label='ReX')
    plt.plot(omega_w, Xt_w.imag, label='ImXt')
    plt.plot(omega_w, Xt_w.real, label='ReXt')
    plt.plot(omega_w, Xh_w.imag, label='ImXh')
    plt.plot(omega_w, Xh_w.real, label='ReXh')
    plt.legend()
    plt.show()
