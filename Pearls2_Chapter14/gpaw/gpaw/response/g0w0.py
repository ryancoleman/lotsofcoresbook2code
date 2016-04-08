# This makes sure that division works as in Python3
from __future__ import division, print_function

import functools
import pickle
from math import pi

import numpy as np
from ase.dft.kpoints import monkhorst_pack
from ase.units import Hartree
from ase.utils import opencew, devnull
from ase.utils.timing import timer

import gpaw.mpi as mpi
from gpaw import debug
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.response.chi0 import Chi0, HilbertTransform
from gpaw.response.pair import PairDensity
from gpaw.response.wstc import WignerSeitzTruncatedCoulomb
from gpaw.wavefunctions.pw import PWDescriptor, count_reciprocal_vectors
from gpaw.xc.exx import EXX, select_kpts
from gpaw.xc.tools import vxc


class G0W0(PairDensity):
    """This class defines the G0W0 calculator. The G0W0 calculator is used
    is used to calculate the quasi particle energies through the G0W0
    approximation for a number of states.

    Note: So far the G0W0 calculation only works for spin-paired systems.

    Parameters:
       calc: str or PAW object
          GPAW calculator object or filename of saved calculator object.
       filename: str
          Base filename of output files.
       kpts: list
          List of indices of the IBZ k-points to calculate the quasi particle
          energies for.
       bands: tuple
          Range of band indices, like (n1, n2+1), to calculate the quasi
          particle energies for. Note that the second band index is not
          included.
       ecut: float
          Plane wave cut-off energy in eV.
       nbands: int
          Number of bands to use in the calculation. If :None: the number will
          be determined from :ecut: to yield a number close to the number of
          plane waves used.
       ppa: bool
          Sets whether the Godby-Needs plasmon-pole approximation for the
          dielectric function should be used.
       E0: float
          Energy (in eV) used for fitting in the plasmon-pole approximation.
       domega0: float
          Minimum frequency step (in eV) used in the generation of the non-
          linear frequency grid.
       omega2: float
          Control parameter for the non-linear frequency grid, equal to the
          frequency where the grid spacing has doubled in size.
       wstc: bool
          Sets whether a Wigner-Seitz truncation should be used for the
          Coloumb potential.
    """
    def __init__(self, calc, filename='gw',
                 kpts=None, bands=None, nbands=None, ppa=False,
                 wstc=False,
                 ecut=150.0, eta=0.1, E0=1.0 * Hartree,
                 domega0=0.025, omega2=10.0,
                 nblocks=1, savew=False,
                 world=mpi.world):

        if world.rank != 0:
            txt = devnull
        else:
            txt = open(filename + '.txt', 'w')

        p = functools.partial(print, file=txt)
        p('  ___  _ _ _ ')
        p(' |   || | | |')
        p(' | | || | | |')
        p(' |__ ||_____|')
        p(' |___|')
        p()

        PairDensity.__init__(self, calc, ecut, world=world, nblocks=nblocks,
                             txt=txt)

        self.filename = filename
        self.savew = savew
        
        ecut /= Hartree
        
        self.ppa = ppa
        self.wstc = wstc
        self.eta = eta / Hartree
        self.E0 = E0 / Hartree
        self.domega0 = domega0 / Hartree
        self.omega2 = omega2 / Hartree

        self.kpts = select_kpts(kpts, self.calc)
                
        if bands is None:
            bands = [0, self.nocc2]
            
        self.bands = bands

        b1, b2 = bands
        self.shape = shape = (self.calc.wfs.nspins, len(self.kpts), b2 - b1)
        self.eps_sin = np.empty(shape)     # KS-eigenvalues
        self.f_sin = np.empty(shape)       # occupation numbers
        self.sigma_sin = np.zeros(shape)   # self-energies
        self.dsigma_sin = np.zeros(shape)  # derivatives of self-energies
        self.vxc_sin = None                # KS XC-contributions
        self.exx_sin = None                # exact exchange contributions
        self.Z_sin = None                  # renormalization factors
        
        if nbands is None:
            nbands = int(self.vol * ecut**1.5 * 2**0.5 / 3 / pi**2)
        self.nbands = nbands

        p()
        p('Quasi particle states:')
        if kpts is None:
            p('All k-points in IBZ')
        else:
            kptstxt = ', '.join(['{0:d}'.format(k) for k in self.kpts])
            p('k-points (IBZ indices): [' + kptstxt + ']')
        p('Band range: ({0:d}, {1:d})'.format(b1, b2))
        p()
        p('Computational parameters:')
        p('Plane wave cut-off: {0:g} eV'.format(self.ecut * Hartree))
        p('Number of bands: {0:d}'.format(self.nbands))
        p('Broadening: {0:g} eV'.format(self.eta * Hartree))

        kd = self.calc.wfs.kd

        self.mysKn1n2 = None  # my (s, K, n1, n2) indices
        self.distribute_k_points_and_bands(b1, b2, kd.ibz2bz_k[self.kpts])
        
        # Find q-vectors and weights in the IBZ:
        assert -1 not in kd.bz2bz_ks
        offset_c = 0.5 * ((kd.N_c + 1) % 2) / kd.N_c
        bzq_qc = monkhorst_pack(kd.N_c) + offset_c
        self.qd = KPointDescriptor(bzq_qc)
        self.qd.set_symmetry(self.calc.atoms, kd.symmetry)
        
        assert self.calc.wfs.nspins == 1
        
    @timer('G0W0')
    def calculate(self, ecuts=None):
        """Starts the G0W0 calculation. Returns a dict with the results with
        the following key/value pairs:

        f: (s, k, n) ndarray
           Occupation numbers
        eps: (s, k, n) ndarray
           Kohn-Sham eigenvalues in eV
        vxc: (s, k, n) ndarray
           Exchange-correlation contributions in eV
        exx: (s, k, n) ndarray
           Exact exchange contributions in eV
        sigma: (s, k, n) ndarray
           Self-energy contributions in eV
        Z: (s, k, n) ndarray
           Renormalization factors
        qp: (s, k, n) ndarray
           Quasi particle energies in eV
        """
        kd = self.calc.wfs.kd

        self.calculate_ks_xc_contribution()
        self.calculate_exact_exchange()

        # Reset calculation
        self.sigma_sin = np.zeros(self.shape)   # self-energies
        self.dsigma_sin = np.zeros(self.shape)  # derivatives of self-energies

        # Get KS eigenvalues and occupation numbers:
        b1, b2 = self.bands
        for i, k in enumerate(self.kpts):
            kpt = self.calc.wfs.kpt_u[k]
            self.eps_sin[0, i] = kpt.eps_n[b1:b2]
            self.f_sin[0, i] = kpt.f_n[b1:b2] / kpt.weight

        # My part of the states we want to calculate QP-energies for:
        mykpts = [self.get_k_point(s, K, n1, n2)
                  for s, K, n1, n2 in self.mysKn1n2]
        
        # Loop over q in the IBZ:
        for pd0, W0, q_c in self.calculate_screened_potential():
            for kpt1 in mykpts:
                K2 = kd.find_k_plus_q(q_c, [kpt1.K])[0]
                kpt2 = self.get_k_point(0, K2, 0, self.nbands, block=True)
                k1 = kd.bz2ibz_k[kpt1.K]
                i = self.kpts.index(k1)
                self.calculate_q(i, kpt1, kpt2, pd0, W0)

        self.world.sum(self.sigma_sin)
        self.world.sum(self.dsigma_sin)
        
        self.Z_sin = 1 / (1 - self.dsigma_sin)
        self.qp_sin = self.eps_sin + self.Z_sin * (self.sigma_sin +
                                                   self.exx_sin -
                                                   self.vxc_sin)
        
        results = {'f': self.f_sin,
                   'eps': self.eps_sin * Hartree,
                   'vxc': self.vxc_sin * Hartree,
                   'exx': self.exx_sin * Hartree,
                   'sigma': self.sigma_sin * Hartree,
                   'Z': self.Z_sin,
                   'qp': self.qp_sin * Hartree}
      
        self.print_results(results)
        
        return results
        
    def calculate_q(self, i, kpt1, kpt2, pd0, W0):
        """Calculates the contribution to the self-energy and its derivative
        for a given set of k-points, kpt1 and kpt2."""
        wfs = self.calc.wfs
        
        N_c = pd0.gd.N_c
        i_cG = self.sign * np.dot(self.U_cc,
                                  np.unravel_index(pd0.Q_qG[0], N_c))

        q_c = wfs.kd.bzk_kc[kpt2.K] - wfs.kd.bzk_kc[kpt1.K]
        q0 = np.allclose(q_c, 0) and not self.wstc

        shift0_c = q_c - self.sign * np.dot(self.U_cc, pd0.kd.bzk_kc[0])
        assert np.allclose(shift0_c.round(), shift0_c)
        shift0_c = shift0_c.round().astype(int)

        shift_c = kpt1.shift_c - kpt2.shift_c - shift0_c
        I_G = np.ravel_multi_index(i_cG + shift_c[:, None], N_c, 'wrap')
        
        G_Gv = pd0.get_reciprocal_vectors()
        pos_av = np.dot(self.spos_ac, pd0.gd.cell_cv)
        M_vv = np.dot(pd0.gd.cell_cv.T,
                      np.dot(self.U_cc.T,
                             np.linalg.inv(pd0.gd.cell_cv).T))
        
        Q_aGii = []
        for a, Q_Gii in enumerate(self.Q_aGii):
            x_G = np.exp(1j * np.dot(G_Gv, (pos_av[a] -
                                            self.sign *
                                            np.dot(M_vv, pos_av[a]))))
            U_ii = self.calc.wfs.setups[a].R_sii[self.s]
            Q_Gii = np.dot(np.dot(U_ii, Q_Gii * x_G[:, None, None]),
                           U_ii.T).transpose(1, 0, 2)
            Q_aGii.append(Q_Gii)

        if debug:
            self.check(i_cG, shift0_c, N_c, q_c, Q_aGii)
                
        if self.ppa:
            calculate_sigma = self.calculate_sigma_ppa
        else:
            calculate_sigma = self.calculate_sigma
            
        for n in range(kpt1.n2 - kpt1.n1):
            ut1cc_R = kpt1.ut_nR[n].conj()
            eps1 = kpt1.eps_n[n]
            C1_aGi = [np.dot(Qa_Gii, P1_ni[n].conj())
                      for Qa_Gii, P1_ni in zip(Q_aGii, kpt1.P_ani)]
            n_mG = self.calculate_pair_densities(ut1cc_R, C1_aGi, kpt2,
                                                 pd0, I_G)
            if self.sign == 1:
                n_mG = n_mG.conj()
                
            if q0:
                n_mG[:, 0] = 0
                m = n + kpt1.n1 - kpt2.n1
                if 0 <= m < len(n_mG):
                    n_mG[m, 0] = 1.0
                    
            f_m = kpt2.f_n
            deps_m = eps1 - kpt2.eps_n
            sigma, dsigma = calculate_sigma(n_mG, deps_m, f_m, W0)
            nn = kpt1.n1 + n - self.bands[0]
            self.sigma_sin[kpt1.s, i, nn] += sigma
            self.dsigma_sin[kpt1.s, i, nn] += dsigma
            
    def check(self, i_cG, shift0_c, N_c, q_c, Q_aGii):
        I0_G = np.ravel_multi_index(i_cG - shift0_c[:, None], N_c, 'wrap')
        qd1 = KPointDescriptor([q_c])
        pd1 = PWDescriptor(self.ecut, self.calc.wfs.gd, complex, qd1)
        G_I = np.empty(N_c.prod(), int)
        G_I[:] = -1
        I1_G = pd1.Q_qG[0]
        G_I[I1_G] = np.arange(len(I0_G))
        G_G = G_I[I0_G]
        assert len(I0_G) == len(I1_G)
        assert (G_G >= 0).all()

        for a, Q_Gii in enumerate(self.initialize_paw_corrections(pd1)):
            e = abs(Q_aGii[a] - Q_Gii[G_G]).max()
            assert e < 1e-12

    @timer('Sigma')
    def calculate_sigma(self, n_mG, deps_m, f_m, C_swGG):
        """Calculates a contribution to the self-energy and its derivative for
        a given (k, k-q)-pair from its corresponding pair-density and
        energy."""
        o_m = abs(deps_m)
        # Add small number to avoid zeros for degenerate states:
        sgn_m = np.sign(deps_m + 1e-15)
        
        # Pick +i*eta or -i*eta:
        s_m = (1 + sgn_m * np.sign(0.5 - f_m)).astype(int) // 2
        
        beta = (2**0.5 - 1) * self.domega0 / self.omega2
        w_m = (o_m / (self.domega0 + beta * o_m)).astype(int)
        o1_m = self.omega_w[w_m]
        o2_m = self.omega_w[w_m + 1]
        
        x = 1.0 / (self.qd.nbzkpts * 2 * pi * self.vol)
        sigma = 0.0
        dsigma = 0.0
        # Performing frequency integration
        for o, o1, o2, sgn, s, w, n_G in zip(o_m, o1_m, o2_m,
                                             sgn_m, s_m, w_m, n_mG):
            C1_GG = C_swGG[s][w]
            C2_GG = C_swGG[s][w + 1]
            p = x * sgn
            myn_G = n_G[self.Ga:self.Gb]
            sigma1 = p * np.dot(np.dot(myn_G, C1_GG), n_G.conj()).imag
            sigma2 = p * np.dot(np.dot(myn_G, C2_GG), n_G.conj()).imag
            sigma += ((o - o1) * sigma2 + (o2 - o) * sigma1) / (o2 - o1)
            dsigma += sgn * (sigma2 - sigma1) / (o2 - o1)
            
        return sigma, dsigma

    def calculate_screened_potential(self):
        """Calculates the screened potential for each q-point in the 1st BZ.
        Since many q-points are related by symmetry, the actual calculation is
        only done for q-points in the IBZ and the rest are obtained by symmetry
        transformations. Results are returned as a generator to that it is not
        necessary to store a huge matrix for each q-point in the memory."""
        # The decorator $timer('W') doesn't work for generators, do we will
        # have to manually start and stop the timer here:
        self.timer.start('W')
        print('Calculating screened Coulomb potential', file=self.fd)
        if self.wstc:
            print('Using Wigner-Seitz truncated Coloumb potential',
                  file=self.fd)
            
        if self.ppa:
            print('Using Godby-Needs plasmon-pole approximation:',
                  file=self.fd)
            print('    Fitting energy: i*E0, E0 = %.3f Hartee' % self.E0,
                  file=self.fd)

            # use small imaginary frequency to avoid dividing by zero:
            frequencies = [1e-10j, 1j * self.E0 * Hartree]
            
            parameters = {'eta': 0,
                          'hilbert': False,
                          'timeordered': False,
                          'frequencies': frequencies}
        else:
            print('Using full frequency integration:', file=self.fd)
            print('  domega0: {0:g}'.format(self.domega0 * Hartree),
                  file=self.fd)
            print('  omega2: {0:g}'.format(self.omega2 * Hartree),
                  file=self.fd)

            parameters = {'eta': self.eta * Hartree,
                          'hilbert': True,
                          'timeordered': True,
                          'domega0': self.domega0 * Hartree,
                          'omega2': self.omega2 * Hartree}
        
        chi0 = Chi0(self.calc,
                    nbands=self.nbands,
                    ecut=self.ecut * Hartree,
                    intraband=False,
                    real_space_derivatives=False,
                    txt=self.filename + '.w.txt',
                    timer=self.timer,
                    keep_occupied_states=True,
                    nblocks=self.blockcomm.size,
                    no_optical_limit=self.wstc,
                    **parameters)

        if self.wstc:
            wstc = WignerSeitzTruncatedCoulomb(
                self.calc.wfs.gd.cell_cv,
                self.calc.wfs.kd.N_c,
                chi0.fd)
        else:
            wstc = None
        
        self.omega_w = chi0.omega_w
        self.omegamax = chi0.omegamax
        
        htp = HilbertTransform(self.omega_w, self.eta, gw=True)
        htm = HilbertTransform(self.omega_w, -self.eta, gw=True)

        # Find maximum size of chi-0 matrices:
        gd = self.calc.wfs.gd
        nGmax = max(count_reciprocal_vectors(self.ecut, gd, q_c)
                    for q_c in self.qd.ibzk_kc)
        nw = len(self.omega_w)
        
        size = self.blockcomm.size
        mynGmax = (nGmax + size - 1) // size
        mynw = (nw + size - 1) // size
        
        # Allocate memory in the beginning and use for all q:
        A1_x = np.empty(nw * mynGmax * nGmax, complex)
        A2_x = np.empty(max(mynw * nGmax, nw * mynGmax) * nGmax, complex)
        
        # Need to pause the timer in between iterations
        self.timer.stop('W')
        for iq, q_c in enumerate(self.qd.ibzk_kc):
            self.timer.start('W')
            if self.savew:
                wfilename = self.filename + '.w.q%d.pckl' % iq
                fd = opencew(wfilename)
            if self.savew and fd is None:
                # Read screened potential from file
                with open(wfilename) as fd:
                    pd, W = pickle.load(fd)
            else:
                # First time calculation
                pd, W = self.calculate_w(chi0, q_c, htp, htm, wstc, A1_x, A2_x)
                if self.savew:
                    pickle.dump((pd, W), fd, pickle.HIGHEST_PROTOCOL)

            self.timer.stop('W')
            # Loop over all k-points in the BZ and find those that are related
            # to the current IBZ k-point by symmetry
            Q1 = self.qd.ibz2bz_k[iq]
            done = set()
            for s, Q2 in enumerate(self.qd.bz2bz_ks[Q1]):
                if Q2 >= 0 and Q2 not in done:
                    s = self.qd.sym_k[Q2]
                    self.s = s
                    self.U_cc = self.qd.symmetry.op_scc[s]
                    time_reversal = self.qd.time_reversal_k[Q2]
                    self.sign = 1 - 2 * time_reversal
                    Q_c = self.qd.bzk_kc[Q2]
                    d_c = self.sign * np.dot(self.U_cc, q_c) - Q_c
                    assert np.allclose(d_c.round(), d_c)
                    yield pd, W, Q_c
                    done.add(Q2)
    
    @timer('WW')
    def calculate_w(self, chi0, q_c, htp, htm, wstc, A1_x, A2_x):
        """Calculates the screened potential for a specified q-point."""
        pd, chi0_wGG = chi0.calculate(q_c, A_x=A1_x)[:2]
        self.Q_aGii = chi0.Q_aGii
        self.Ga = chi0.Ga
        self.Gb = chi0.Gb
        
        if self.blockcomm.size > 1:
            A1_x = chi0_wGG.ravel()
            chi0_wGG = chi0.redistribute(chi0_wGG, A2_x)
            
        if self.wstc:
            iG_G = (wstc.get_potential(pd) / (4 * pi))**0.5
            if np.allclose(q_c, 0):
                chi0_wGG[:, 0] = 0.0
                chi0_wGG[:, :, 0] = 0.0
                G0inv = 0.0
                G20inv = 0.0
            else:
                G0inv = None
                G20inv = None
        else:
            if np.allclose(q_c, 0):
                dq3 = (2 * pi)**3 / (self.qd.nbzkpts * self.vol)
                qc = (dq3 / 4 / pi * 3)**(1 / 3)
                G0inv = 2 * pi * qc**2 / dq3
                G20inv = 4 * pi * qc / dq3
                G_G = pd.G2_qG[0]**0.5
                G_G[0] = 1
                iG_G = 1 / G_G
            else:
                iG_G = pd.G2_qG[0]**-0.5
                G0inv = None
                G20inv = None

        delta_GG = np.eye(len(iG_G))

        if self.ppa:
            return pd, self.ppa_w(chi0_wGG, iG_G, delta_GG, G0inv, G20inv, q_c)
            
        self.timer.start('Dyson eq.')
        # Calculate W and store it in chi0_wGG ndarray:
        for chi0_GG in chi0_wGG:
            e_GG = (delta_GG -
                    4 * pi * chi0_GG * iG_G * iG_G[:, np.newaxis])
            W_GG = chi0_GG
            W_GG[:] = 4 * pi * (np.linalg.inv(e_GG) -
                                delta_GG) * iG_G * iG_G[:, np.newaxis]
            if np.allclose(q_c, 0):
                W_GG[0, 0] *= G20inv
                W_GG[1:, 0] *= G0inv
                W_GG[0, 1:] *= G0inv
                
        if self.blockcomm.size > 1:
            Wm_wGG = chi0.redistribute(chi0_wGG, A1_x)
        else:
            Wm_wGG = chi0_wGG
            
        Wp_wGG = A2_x[:Wm_wGG.size].reshape(Wm_wGG.shape)
        Wp_wGG[:] = Wm_wGG

        with self.timer('Hilbert transform'):
            htp(Wp_wGG)
            htm(Wm_wGG)
        self.timer.stop('Dyson eq.')
        
        return pd, [Wp_wGG, Wm_wGG]

    @timer('Kohn-Sham XC-contribution')
    def calculate_ks_xc_contribution(self):
        name = self.filename + '.vxc.npy'
        fd = opencew(name)
        if fd is None:
            print('Reading Kohn-Sham XC contribution from file:', name,
                  file=self.fd)
            with open(name) as fd:
                self.vxc_sin = np.load(fd)
            assert self.vxc_sin.shape == self.shape, self.vxc_sin.shape
            return
            
        print('Calculating Kohn-Sham XC contribution', file=self.fd)
        if self.reader is not None:
            self.calc.wfs.read_projections(self.reader)
        vxc_skn = vxc(self.calc, self.calc.hamiltonian.xc) / Hartree
        n1, n2 = self.bands
        self.vxc_sin = vxc_skn[:, self.kpts, n1:n2]
        np.save(fd, self.vxc_sin)
        
    @timer('EXX')
    def calculate_exact_exchange(self):
        name = self.filename + '.exx.npy'
        fd = opencew(name)
        if fd is None:
            print('Reading EXX contribution from file:', name, file=self.fd)
            with open(name) as fd:
                self.exx_sin = np.load(fd)
            assert self.exx_sin.shape == self.shape, self.exx_sin.shape
            return
            
        print('Calculating EXX contribution', file=self.fd)
        exx = EXX(self.calc, kpts=self.kpts, bands=self.bands,
                  txt=self.filename + '.exx.txt', timer=self.timer)
        exx.calculate()
        self.exx_sin = exx.get_eigenvalue_contributions() / Hartree
        np.save(fd, self.exx_sin)

    def print_results(self, results):
        description = ['f:     Occupation numbers',
                       'eps:   KS-eigenvalues [eV]',
                       'vxc:   KS vxc [eV]',
                       'exx:   Exact exchange [eV]',
                       'sigma: Self-energies [eV]',
                       'Z:     Renormalization factors',
                       'qp:    QP-energies [eV]']

        print('\nResults:', file=self.fd)
        for line in description:
            print(line, file=self.fd)
            
        b1, b2 = self.bands
        names = [line.split(':', 1)[0] for line in description]
        ibzk_kc = self.calc.wfs.kd.ibzk_kc
        for s in range(self.calc.wfs.nspins):
            for i, ik in enumerate(self.kpts):
                print('\nk-point ' +
                      '{0} ({1}): ({2:.3f}, {3:.3f}, {4:.3f})'.format(
                          i, ik, *ibzk_kc[ik]), file=self.fd)
                print('band' +
                      ''.join('{0:>8}'.format(name) for name in names),
                      file=self.fd)
                for n in range(b2 - b1):
                    print('{0:4}'.format(n + b1) +
                          ''.join('{0:8.3f}'.format(results[name][s, i, n])
                                  for name in names),
                          file=self.fd)

        self.timer.write(self.fd)

    @timer('PPA')
    def ppa_w(self, chi0_wGG, iG_G, delta_GG, G0inv, G20inv, q_c):
        einv_wGG = []
        for chi0_GG in chi0_wGG:
            e_GG = (delta_GG -
                    4 * pi * chi0_GG * iG_G * iG_G[:, np.newaxis])
            einv_wGG.append(np.linalg.inv(e_GG) - delta_GG)

        if self.wstc and np.allclose(q_c, 0):
            einv_wGG[0][0] = 42
            einv_wGG[0][:, 0] = 42
        omegat_GG = self.E0 * np.sqrt(einv_wGG[1] /
                                      (einv_wGG[0] - einv_wGG[1]))
        R_GG = -0.5 * omegat_GG * einv_wGG[0]
        W_GG = 4 * pi**2 * R_GG * iG_G * iG_G[:, np.newaxis]
        if np.allclose(q_c, 0):
            W_GG[0, 0] *= G20inv
            W_GG[1:, 0] *= G0inv
            W_GG[0, 1:] *= G0inv

        return [W_GG, omegat_GG]
        
    @timer('PPA-Sigma')
    def calculate_sigma_ppa(self, n_mG, deps_m, f_m, W):
        W_GG, omegat_GG = W

        sigma = 0.0
        dsigma = 0.0

        # init variables (is this necessary?)
        nG = n_mG.shape[1]
        deps_GG = np.empty((nG, nG))
        sign_GG = np.empty((nG, nG))
        x1_GG = np.empty((nG, nG))
        x2_GG = np.empty((nG, nG))
        x3_GG = np.empty((nG, nG))
        x4_GG = np.empty((nG, nG))
        x_GG = np.empty((nG, nG))
        dx_GG = np.empty((nG, nG))
        nW_G = np.empty(nG)
        for m in range(np.shape(n_mG)[0]):
            deps_GG = deps_m[m]
            sign_GG = 2 * f_m[m] - 1
            x1_GG = 1 / (deps_GG + omegat_GG - 1j * self.eta)
            x2_GG = 1 / (deps_GG - omegat_GG + 1j * self.eta)
            x3_GG = 1 / (deps_GG + omegat_GG - 1j * self.eta * sign_GG)
            x4_GG = 1 / (deps_GG - omegat_GG - 1j * self.eta * sign_GG)
            x_GG = W_GG * (sign_GG * (x1_GG - x2_GG) + x3_GG + x4_GG)
            dx_GG = W_GG * (sign_GG * (x1_GG**2 - x2_GG**2) +
                            x3_GG**2 + x4_GG**2)
            nW_G = np.dot(n_mG[m], x_GG)
            sigma += np.vdot(n_mG[m], nW_G).real
            nW_G = np.dot(n_mG[m], dx_GG)
            dsigma -= np.vdot(n_mG[m], nW_G).real
        
        x = 1 / (self.qd.nbzkpts * 2 * pi * self.vol)
        return x * sigma, x * dsigma
