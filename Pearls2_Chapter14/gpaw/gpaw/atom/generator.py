# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.
from __future__ import print_function
from math import pi, sqrt

import numpy as np
from numpy.linalg import solve, inv
from ase.data import atomic_names

from gpaw.setup_data import SetupData
from gpaw.atom.configurations import configurations
from gpaw.version import version
from gpaw.atom.all_electron import AllElectron, shoot
from gpaw.utilities.lapack import general_diagonalize
from gpaw.utilities import hartree
from gpaw.xc.hybrid import constructX, atomic_exact_exchange
from gpaw.atom.filter import Filter


class Generator(AllElectron):
    def __init__(self, symbol, xcname='LDA', scalarrel=False, corehole=None,
                 configuration=None,
                 nofiles=True, txt='-', gpernode=150,
                 orbital_free=False, tf_coeff=1.):
        AllElectron.__init__(self, symbol, xcname, scalarrel, corehole,
                             configuration, nofiles, txt, gpernode,
                             orbital_free, tf_coeff)

    def run(self, core='', rcut=1.0, extra=None,
            logderiv=False, vbar=None, exx=False, name=None,
            normconserving='', filter=(0.4, 1.75), rcutcomp=None,
            write_xml=True, use_restart_file=True,
            empty_states=''):

        self.name = name

        self.core = core
        if type(rcut) is float:
            rcut_l = [rcut]
        else:
            rcut_l = rcut
        rcutmax = max(rcut_l)
        rcutmin = min(rcut_l)
        self.rcut_l = rcut_l

        if rcutcomp is None:
            rcutcomp = rcutmin
        self.rcutcomp = rcutcomp

        hfilter, xfilter = filter

        Z = self.Z

        n_j = self.n_j
        l_j = self.l_j
        f_j = self.f_j
        e_j = self.e_j

        if vbar is None:
            vbar = ('poly', rcutmin * 0.9)
        vbar_type, rcutvbar = vbar

        normconserving_l = [x in normconserving for x in 'spdf']

        # Parse core string:
        j = 0
        if core.startswith('['):
            a, core = core.split(']')
            core_symbol = a[1:]
            j = len(configurations[core_symbol][1])

        while core != '':
            assert n_j[j] == int(core[0])
            assert l_j[j] == 'spdf'.find(core[1])
            if j != self.jcorehole:
                assert f_j[j] == 2 * (2 * l_j[j] + 1)
            j += 1
            core = core[2:]

        njcore = j
        self.njcore = njcore

        while empty_states != '':
            n = int(empty_states[0])
            l = 'spdf'.find(empty_states[1])
            assert n == 1 + l + l_j.count(l)
            n_j.append(n)
            l_j.append(l)
            f_j.append(0.0)
            e_j.append(-0.01)
            empty_states = empty_states[2:]
            
        if 2 in l_j[njcore:]:
            # We have a bound valence d-state.  Add bound s- and
            # p-states if not already there:
            for l in [0, 1]:
                if l not in l_j[njcore:]:
                    n_j.append(1 + l + l_j.count(l))
                    l_j.append(l)
                    f_j.append(0.0)
                    e_j.append(-0.01)

        if l_j[njcore:] == [0] and Z > 2:
            # We have only a bound valence s-state and we are not
            # hydrogen and not helium.  Add bound p-state:
            n_j.append(n_j[njcore])
            l_j.append(1)
            f_j.append(0.0)
            e_j.append(-0.01)

        nj = len(n_j)

        self.Nv = sum(f_j[njcore:])
        self.Nc = sum(f_j[:njcore])

        lmaxocc = max(l_j[njcore:])
        lmax = max(l_j[njcore:])

        #Parameters for orbital_free
        if self.orbital_free:
            self.n_j = [1]
            self.l_j = [0]
            self.f_j = [self.Z]
            self.e_j = [self.e_j[0]]
            n_j = self.n_j
            l_j = self.l_j
            f_j = self.f_j
            e_j = self.e_j
            nj = len(n_j)
            lmax = 0
            lmaxocc = 0

        # Do all-electron calculation:
        AllElectron.run(self, use_restart_file)

        # Highest occupied atomic orbital:
        self.emax = max(e_j)

        N = self.N
        r = self.r
        dr = self.dr
        d2gdr2 = self.d2gdr2
        beta = self.beta

        dv = r**2 * dr

        t = self.text
        t()
        t('Generating PAW setup')
        if core != '':
            t('Frozen core:', core)

        # So far - no ghost-states:
        self.ghost = False

        # Calculate the kinetic energy of the core states:
        Ekincore = 0.0
        j = 0
        for f, e, u in zip(f_j[:njcore], e_j[:njcore], self.u_j[:njcore]):
            u = np.where(abs(u) < 1e-160, 0, u)  # XXX Numeric!
            k = e - np.sum((u**2 * self.vr * dr)[1:] / r[1:])
            Ekincore += f * k
            if j == self.jcorehole:
                self.Ekincorehole = k
            j += 1

        # Calculate core density:
        if njcore == 0:
            nc = np.zeros(N)
        else:
            uc_j = self.u_j[:njcore]
            uc_j = np.where(abs(uc_j) < 1e-160, 0, uc_j)  # XXX Numeric!
            nc = np.dot(f_j[:njcore], uc_j**2) / (4 * pi)
            nc[1:] /= r[1:]**2
            nc[0] = nc[1]

        self.nc = nc

        # Calculate core kinetic energy density
        if njcore == 0:
            tauc = np.zeros(N)
        else:
            tauc = self.radial_kinetic_energy_density(f_j[:njcore],
                                                      l_j[:njcore],
                                                      self.u_j[:njcore])
            t('Kinetic energy of the core from tauc =',
              np.dot(tauc * r * r, dr) * 4 * pi)

        # Order valence states with respect to angular momentum
        # quantum number:
        self.n_ln = n_ln = []
        self.f_ln = f_ln = []
        self.e_ln = e_ln = []
        for l in range(lmax + 1):
            n_n = []
            f_n = []
            e_n = []
            for j in range(njcore, nj):
                if l_j[j] == l:
                    n_n.append(n_j[j])
                    f_n.append(f_j[j])
                    e_n.append(e_j[j])
            n_ln.append(n_n)
            f_ln.append(f_n)
            e_ln.append(e_n)

        # Add extra projectors:
        if extra is not None:
            if len(extra) == 0:
                lmaxextra = 0
            else:
                lmaxextra = max(extra.keys())
            if lmaxextra > lmax:
                for l in range(lmax, lmaxextra):
                    n_ln.append([])
                    f_ln.append([])
                    e_ln.append([])
                lmax = lmaxextra
            for l in extra:
                nn = -1
                for e in extra[l]:
                    n_ln[l].append(nn)
                    f_ln[l].append(0.0)
                    e_ln[l].append(e)
                    nn -= 1
        else:
            # Automatic:

            # Make sure we have two projectors for each occupied channel:
            for l in range(lmaxocc + 1):
                if len(n_ln[l]) < 2 and not normconserving_l[l]:
                    # Only one - add one more:
                    assert len(n_ln[l]) == 1
                    n_ln[l].append(-1)
                    f_ln[l].append(0.0)
                    e_ln[l].append(1.0 + e_ln[l][0])

            if lmaxocc < 2 and lmaxocc == lmax:
                # Add extra projector for l = lmax + 1:
                n_ln.append([-1])
                f_ln.append([0.0])
                e_ln.append([0.0])
                lmax += 1

        self.lmax = lmax

        rcut_l.extend([rcutmin] * (lmax + 1 - len(rcut_l)))

        t('Cutoffs:')
        for rc, s in zip(rcut_l, 'spdf'):
            t('rc(%s)=%.3f' % (s, rc))
        t('rc(vbar)=%.3f' % rcutvbar)
        t('rc(comp)=%.3f' % rcutcomp)
        t('rc(nct)=%.3f' % rcutmax)
        t()
        t('Kinetic energy of the core states: %.6f' % Ekincore)

        # Allocate arrays:
        self.u_ln = u_ln = []  # phi * r
        self.s_ln = s_ln = []  # phi-tilde * r
        self.q_ln = q_ln = []  # p-tilde * r
        for l in range(lmax + 1):
            nn = len(n_ln[l])
            u_ln.append(np.zeros((nn, N)))
            s_ln.append(np.zeros((nn, N)))
            q_ln.append(np.zeros((nn, N)))

        # Fill in all-electron wave functions:
        for l in range(lmax + 1):
            # Collect all-electron wave functions:
            u_n = [self.u_j[j] for j in range(njcore, nj) if l_j[j] == l]
            for n, u in enumerate(u_n):
                u_ln[l][n] = u

        # Grid-index corresponding to rcut:
        gcut_l = [1 + int(rc * N / (rc + beta)) for rc in rcut_l]

        rcutfilter = xfilter * rcutmax
        self.rcutfilter = rcutfilter
        gcutfilter = 1 + int(rcutfilter * N / (rcutfilter + beta))
        gcutmax = 1 + int(rcutmax * N / (rcutmax + beta))

        # Outward integration of unbound states stops at 3 * rcut:
        gmax = int(3 * rcutmax * N / (3 * rcutmax + beta))
        assert gmax > gcutfilter

        # Calculate unbound extra states:
        c2 = -(r / dr)**2
        c10 = -d2gdr2 * r**2
        for l, (n_n, e_n, u_n) in enumerate(zip(n_ln, e_ln, u_ln)):
            for n, e, u in zip(n_n, e_n, u_n):
                if n < 0:
                    u[:] = 0.0
                    shoot(u, l, self.vr, e, self.r2dvdr, r, dr, c10, c2,
                          self.scalarrel, gmax=gmax)
                    u *= 1.0 / u[gcut_l[l]]

        charge = Z - self.Nv - self.Nc
        t('Charge: %.1f' % charge)
        t('Core electrons: %.1f' % self.Nc)
        t('Valence electrons: %.1f' % self.Nv)

        # Construct smooth wave functions:
        coefs = []
        for l, (u_n, s_n) in enumerate(zip(u_ln, s_ln)):
            nodeless = True
            gc = gcut_l[l]
            for u, s in zip(u_n, s_n):
                s[:] = u
                if normconserving_l[l]:
                    A = np.zeros((5, 5))
                    A[:4, 0] = 1.0
                    A[:4, 1] = r[gc - 2:gc + 2]**2
                    A[:4, 2] = A[:4, 1]**2
                    A[:4, 3] = A[:4, 1] * A[:4, 2]
                    A[:4, 4] = A[:4, 2]**2
                    A[4, 4] = 1.0
                    a = u[gc - 2:gc + 3] / r[gc - 2:gc + 3]**(l + 1)
                    a = np.log(a)
                    def f(x):
                        a[4] = x
                        b = solve(A, a)
                        r1 = r[:gc]
                        r2 = r1**2
                        rl1 = r1**(l + 1)
                        y = b[0] + r2 * (b[1] + r2 * (b[2] + r2 * (b[3] + r2
                                                                   * b[4])))
                        y = np.exp(y)
                        s[:gc] = rl1 * y
                        return np.dot(s**2, dr) - 1
                    x1 = 0.0
                    x2 = 0.001
                    f1 = f(x1)
                    f2 = f(x2)
                    while abs(f1) > 1e-6:
                        x0 = (x1 / f1 - x2 / f2) / (1 / f1 - 1 / f2)
                        f0 = f(x0)
                        if abs(f1) < abs(f2):
                            x2, f2 = x1, f1
                        x1, f1 = x0, f0

                else:
                    A = np.ones((4, 4))
                    A[:, 0] = 1.0
                    A[:, 1] = r[gc - 2:gc + 2]**2
                    A[:, 2] = A[:, 1]**2
                    A[:, 3] = A[:, 1] * A[:, 2]
                    a = u[gc - 2:gc + 2] / r[gc - 2:gc + 2]**(l + 1)
                    if 0:#l < 2 and nodeless:
                        a = np.log(a)
                    a = solve(A, a)
                    r1 = r[:gc]
                    r2 = r1**2
                    rl1 = r1**(l + 1)
                    y = a[0] + r2 * (a[1] + r2 * (a[2] + r2 * (a[3])))
                    if 0:#l < 2 and nodeless:
                        y = np.exp(y)
                    s[:gc] = rl1 * y

                coefs.append(a)
                if nodeless:
                    if not np.alltrue(s[1:gc] > 0.0):
                        raise RuntimeError(
                            'Error: The %d%s pseudo wave has a node!' %
                            (n_ln[l][0], 'spdf'[l]))
                    # Only the first state for each l must be nodeless:
                    nodeless = False

        # Calculate pseudo core density:
        gcutnc = 1 + int(rcutmax * N / (rcutmax + beta))
        self.nct = nct = nc.copy()
        A = np.ones((4, 4))
        A[0] = 1.0
        A[1] = r[gcutnc - 2:gcutnc + 2]**2
        A[2] = A[1]**2
        A[3] = A[1] * A[2]
        a = nc[gcutnc - 2:gcutnc + 2]
        a = solve(np.transpose(A), a)
        r2 = r[:gcutnc]**2
        nct[:gcutnc] = a[0] + r2 * (a[1] + r2 * (a[2] + r2 * a[3]))
        t('Pseudo-core charge: %.6f' % (4 * pi * np.dot(nct, dv)))

        # ... and the pseudo core kinetic energy density:
        tauct = tauc.copy()
        a = tauc[gcutnc - 2:gcutnc + 2]
        a = solve(np.transpose(A), a)
        tauct[:gcutnc] = a[0] + r2 * (a[1] + r2 * (a[2] + r2 * a[3]))

        # ... and the soft valence density:
        nt = np.zeros(N)
        for f_n, s_n in zip(f_ln, s_ln):
            nt += np.dot(f_n, s_n**2) / (4 * pi)
        nt[1:] /= r[1:]**2
        nt[0] = nt[1]
        nt += nct
        self.nt = nt

        # Calculate the shape function:
        x = r / rcutcomp
        gaussian = np.zeros(N)
        self.gamma = gamma = 10.0
        gaussian[:gmax] = np.exp(-gamma * x[:gmax]**2)
        gt = 4 * (gamma / rcutcomp**2)**1.5 / sqrt(pi) * gaussian
        t('Shape function alpha=%.3f' % (gamma / rcutcomp**2))
        norm = np.dot(gt, dv)
        #print norm, norm-1
        assert abs(norm - 1) < 1e-2
        gt /= norm

        # Calculate smooth charge density:
        Nt = np.dot(nt, dv)
        rhot = nt - (Nt + charge / (4 * pi)) * gt
        t('Pseudo-electron charge', 4 * pi * Nt)

        vHt = np.zeros(N)
        hartree(0, rhot * r * dr, r, vHt)
        vHt[1:] /= r[1:]
        vHt[0] = vHt[1]

        vXCt = np.zeros(N)

        extra_xc_data = {}

        if self.xc.type != 'GLLB':
            self.xc.calculate_spherical(self.rgd,
                                        nt.reshape((1, -1)),
                                        vXCt.reshape((1, -1)))
        else:
            self.xc.get_smooth_xc_potential_and_energy_1d(vXCt)

            # Calculate extra-stuff for non-local functionals
            self.xc.get_extra_setup_data(extra_xc_data)

        vt = vHt + vXCt

        # Construct zero potential:
        gc = 1 + int(rcutvbar * N / (rcutvbar + beta))
        if vbar_type == 'f':
            assert lmax == 2
            uf = np.zeros(N)
            l = 3

            # Solve for all-electron f-state:
            eps = 0.0
            shoot(uf, l, self.vr, eps, self.r2dvdr, r, dr, c10, c2,
                  self.scalarrel, gmax=gmax)
            uf *= 1.0 / uf[gc]

            # Fit smooth pseudo f-state polynomium:
            A = np.ones((4, 4))
            A[:, 0] = 1.0
            A[:, 1] = r[gc - 2:gc + 2]**2
            A[:, 2] = A[:, 1]**2
            A[:, 3] = A[:, 1] * A[:, 2]
            a = uf[gc - 2:gc + 2] / r[gc - 2:gc + 2]**(l + 1)
            a0, a1, a2, a3 = solve(A, a)
            r1 = r[:gc]
            r2 = r1**2
            rl1 = r1**(l + 1)
            y = a0 + r2 * (a1 + r2 * (a2 + r2 * a3))
            sf = uf.copy()
            sf[:gc] = rl1 * y

            # From 0 to gc, use analytic formula for kinetic energy operator:
            r4 = r2**2
            r6 = r4 * r2
            enumerator = (a0 * l * (l + 1) +
                          a1 * (l + 2) * (l + 3) * r2 +
                          a2 * (l + 4) * (l + 5) * r4 +
                          a3 * (l + 6) * (l + 7) * r6)
            denominator = a0 + a1 * r2 + a2 * r4 + a3 * r6
            ekin_over_phit = - 0.5 * (enumerator / denominator - l * (l + 1))
            ekin_over_phit[1:] /= r2[1:]

            vbar = eps - vt
            vbar[:gc] -= ekin_over_phit
            vbar[0] = vbar[1] # Actually we can collect the terms into
            # a single fraction without poles, so as to avoid doing this,
            # but this is good enough

            # From gc to gmax, use finite-difference formula for kinetic
            # energy operator:
            vbar[gc:gmax] -= self.kin(l, sf)[gc:gmax] / sf[gc:gmax]
            vbar[gmax:] = 0.0
        else:
            assert vbar_type == 'poly'
            A = np.ones((2, 2))
            A[0] = 1.0
            A[1] = r[gc - 1:gc + 1]**2
            a = vt[gc - 1:gc + 1]
            a = solve(np.transpose(A), a)
            r2 = r**2
            vbar = a[0] + r2 * a[1] - vt
            vbar[gc:] = 0.0

        vt += vbar

        # Construct projector functions:
        for l, (e_n, s_n, q_n) in enumerate(zip(e_ln, s_ln, q_ln)):
            for e, s, q in zip(e_n, s_n, q_n):
                q[:] = self.kin(l, s) + (vt - e) * s
                q[gcutmax:] = 0.0

        filter = Filter(r, dr, gcutfilter, hfilter).filter

        vbar = filter(vbar * r)

        # Calculate matrix elements:
        self.dK_lnn = dK_lnn = []
        self.dH_lnn = dH_lnn = []
        self.dO_lnn = dO_lnn = []

        for l, (e_n, u_n, s_n, q_n) in enumerate(zip(e_ln, u_ln,
                                                     s_ln, q_ln)):
            A_nn = np.inner(s_n, q_n * dr)
            # Do a LU decomposition of A:
            nn = len(e_n)
            L_nn = np.identity(nn, float)
            U_nn = A_nn.copy()
            
            # Keep all bound states normalized
            if sum([n > 0 for n in n_ln[l]]) <= 1:
                for i in range(nn):
                    for j in range(i + 1, nn):
                        L_nn[j, i] = 1.0 * U_nn[j, i] / U_nn[i, i]
                        U_nn[j, :] -= U_nn[i, :] * L_nn[j, i]

            dO_nn = (np.inner(u_n, u_n * dr) -
                     np.inner(s_n, s_n * dr))

            e_nn = np.zeros((nn, nn))
            e_nn.ravel()[::nn + 1] = e_n
            dH_nn = np.dot(dO_nn, e_nn) - A_nn

            q_n[:] = np.dot(inv(np.transpose(U_nn)), q_n)
            s_n[:] = np.dot(inv(L_nn), s_n)
            u_n[:] = np.dot(inv(L_nn), u_n)

            dO_nn = np.dot(np.dot(inv(L_nn), dO_nn),
                            inv(np.transpose(L_nn)))
            dH_nn = np.dot(np.dot(inv(L_nn), dH_nn),
                            inv(np.transpose(L_nn)))

            ku_n = [self.kin(l, u, e) for u, e in zip(u_n, e_n)]
            ks_n = [self.kin(l, s) for s in s_n]
            dK_nn = 0.5 * (np.inner(u_n, ku_n * dr) -
                           np.inner(s_n, ks_n * dr))
            dK_nn += np.transpose(dK_nn).copy()

            dK_lnn.append(dK_nn)
            dO_lnn.append(dO_nn)
            dH_lnn.append(dH_nn)

            for n, q in enumerate(q_n):
                q[:] = filter(q, l) * r**(l + 1)

            A_nn = np.inner(s_n, q_n * dr)
            q_n[:] = np.dot(inv(np.transpose(A_nn)), q_n)

        self.vt = vt
        self.vbar = vbar

        t('state    eigenvalue         norm')
        t('--------------------------------')
        for l, (n_n, f_n, e_n) in enumerate(zip(n_ln, f_ln, e_ln)):
            for n in range(len(e_n)):
                if n_n[n] > 0:
                    f = '(%d)' % f_n[n]
                    t('%d%s%-4s: %12.6f %12.6f' % (
                        n_n[n], 'spdf'[l], f, e_n[n],
                        np.dot(s_ln[l][n]**2, dr)))
                else:
                    t('*%s    : %12.6f' % ('spdf'[l], e_n[n]))
        t('--------------------------------')

        self.logd = {}
        if logderiv:
            ni = 300
            self.elog = np.linspace(-5.0, 1.0, ni)
            # Calculate logarithmic derivatives:
            gld = gcutmax + 10
            self.rlog = r[gld]
            assert gld < gmax
            t('Calculating logarithmic derivatives at r=%.3f' % r[gld])
            t('(skip with [Ctrl-C])')

            try:
                u = np.zeros(N)
                for l in range(4):
                    self.logd[l] = (np.empty(ni), np.empty(ni))
                    if l <= lmax:
                        dO_nn = dO_lnn[l]
                        dH_nn = dH_lnn[l]
                        q_n = q_ln[l]

                    fae = open(self.symbol + '.ae.ld.' + 'spdf'[l], 'w')
                    fps = open(self.symbol + '.ps.ld.' + 'spdf'[l], 'w')

                    for i, e in enumerate(self.elog):
                        # All-electron logarithmic derivative:
                        u[:] = 0.0
                        shoot(u, l, self.vr, e, self.r2dvdr, r, dr, c10, c2,
                              self.scalarrel, gmax=gld)
                        dudr = 0.5 * (u[gld + 1] - u[gld - 1]) / dr[gld]
                        ld = dudr / u[gld] - 1.0 / r[gld]
                        print(e, ld, file=fae)
                        self.logd[l][0][i] = ld

                        # PAW logarithmic derivative:
                        s = self.integrate(l, vt, e, gld)
                        if l <= lmax:
                            A_nn = dH_nn - e * dO_nn
                            s_n = [self.integrate(l, vt, e, gld, q)
                                   for q in q_n]
                            B_nn = np.inner(q_n, s_n * dr)
                            a_n = np.dot(q_n, s * dr)

                            B_nn = np.dot(A_nn, B_nn)
                            B_nn.ravel()[::len(a_n) + 1] += 1.0
                            c_n = solve(B_nn, np.dot(A_nn, a_n))
                            s -= np.dot(c_n, s_n)

                        dsdr = 0.5 * (s[gld + 1] - s[gld - 1]) / dr[gld]
                        ld = dsdr / s[gld] - 1.0 / r[gld]
                        print(e, ld, file=fps)
                        self.logd[l][1][i] = ld

            except KeyboardInterrupt:
                pass

        self.write(nc,'nc')
        self.write(nt, 'nt')
        self.write(nct, 'nct')
        self.write(vbar, 'vbar')
        self.write(vt, 'vt')
        self.write(tauc, 'tauc')
        self.write(tauct, 'tauct')

        for l, (n_n, f_n, u_n, s_n, q_n) in enumerate(zip(n_ln, f_ln,
                                                          u_ln, s_ln, q_ln)):
            for n, f, u, s, q in zip(n_n, f_n, u_n, s_n, q_n):
                if n < 0:
                    self.write(u, 'ae', n=n, l=l)
                self.write(s, 'ps', n=n, l=l)
                self.write(q, 'proj', n=n, l=l)

        # Test for ghost states:
        for h in [0.05]:
            self.diagonalize(h)

        self.vn_j = vn_j = []
        self.vl_j = vl_j = []
        self.vf_j = vf_j = []
        self.ve_j = ve_j = []
        self.vu_j = vu_j = []
        self.vs_j = vs_j = []
        self.vq_j = vq_j = []
        j_ln = [[0 for f in f_n] for f_n in f_ln]
        j = 0
        for l, n_n in enumerate(n_ln):
            for n, nn in enumerate(n_n):
                if nn > 0:
                    vf_j.append(f_ln[l][n])
                    vn_j.append(nn)
                    vl_j.append(l)
                    ve_j.append(e_ln[l][n])
                    vu_j.append(u_ln[l][n])
                    vs_j.append(s_ln[l][n])
                    vq_j.append(q_ln[l][n])
                    j_ln[l][n] = j
                    j += 1
        for l, n_n in enumerate(n_ln):
            for n, nn in enumerate(n_n):
                if nn < 0:
                    vf_j.append(0)
                    vn_j.append(nn)
                    vl_j.append(l)
                    ve_j.append(e_ln[l][n])
                    vu_j.append(u_ln[l][n])
                    vs_j.append(s_ln[l][n])
                    vq_j.append(q_ln[l][n])
                    j_ln[l][n] = j
                    j += 1
        nj = j

        self.dK_jj = np.zeros((nj, nj))
        for l, j_n in enumerate(j_ln):
            for n1, j1 in enumerate(j_n):
                for n2, j2 in enumerate(j_n):
                    self.dK_jj[j1, j2] = self.dK_lnn[l][n1, n2]

        if exx:
            X_p = constructX(self)
            ExxC = atomic_exact_exchange(self, 'core-core')
        else:
            X_p = None
            ExxC = None

        sqrt4pi = sqrt(4 * pi)
        setup = SetupData(self.symbol, self.xc.name, self.name,
                          readxml=False)

        def divide_by_r(x_g, l):
            r = self.r
            #for x_g, l in zip(x_jg, l_j):
            p = x_g.copy()
            p[1:] /= self.r[1:]
            # XXXXX go to higher order!!!!!
            if l == 0:#l_j[self.jcorehole] == 0:
                p[0] = (p[2] +
                        (p[1] - p[2]) * (r[0] - r[2]) / (r[1] - r[2]))
            return p

        def divide_all_by_r(x_jg):
            return [divide_by_r(x_g, l) for x_g, l in zip(x_jg, vl_j)]

        setup.l_j = vl_j
        setup.n_j = vn_j
        setup.f_j = vf_j
        setup.eps_j = ve_j
        setup.rcut_j = [rcut_l[l] for l in vl_j]

        setup.nc_g = nc * sqrt4pi
        setup.nct_g = nct * sqrt4pi
        setup.nvt_g = (nt - nct) * sqrt4pi
        setup.e_kinetic_core = Ekincore
        setup.vbar_g = vbar * sqrt4pi
        setup.tauc_g = tauc * sqrt4pi
        setup.tauct_g = tauct * sqrt4pi
        setup.extra_xc_data = extra_xc_data
        setup.Z = Z
        setup.Nc = self.Nc
        setup.Nv = self.Nv
        setup.e_kinetic = self.Ekin
        setup.e_xc = self.Exc
        setup.e_electrostatic = self.Epot
        setup.e_total = self.Epot + self.Exc + self.Ekin

        setup.rgd = self.rgd
        
        setup.rcgauss = self.rcutcomp / sqrt(self.gamma)
        setup.e_kin_jj = self.dK_jj
        setup.ExxC = ExxC
        setup.phi_jg = divide_all_by_r(vu_j)
        setup.phit_jg = divide_all_by_r(vs_j)
        setup.pt_jg = divide_all_by_r(vq_j)
        setup.X_p = X_p

        if self.jcorehole is not None:
            setup.has_corehole = True
            setup.lcorehole = l_j[self.jcorehole] # l_j or vl_j ????? XXX
            setup.ncorehole = n_j[self.jcorehole]
            setup.phicorehole_g = divide_by_r(self.u_j[self.jcorehole],
                                                  setup.lcorehole)
            setup.core_hole_e = self.e_j[self.jcorehole]
            setup.core_hole_e_kin = self.Ekincorehole
            setup.fcorehole = self.fcorehole

        if self.ghost and not self.orbital_free: #In orbital_free we are not interested in ghosts
            raise RuntimeError('Ghost!')

        if self.scalarrel:
            reltype = 'scalar-relativistic'
        else:
            reltype = 'non-relativistic'

        attrs = [('type', reltype), ('name', 'gpaw-%s' % version)]
        data = 'Frozen core: '+ (self.core or 'none')

        setup.generatorattrs = attrs
        setup.generatordata  = data
        setup.orbital_free = self.orbital_free
        
        self.id_j = []
        for l, n in zip(vl_j, vn_j):
            if n > 0:
                id = '%s-%d%s' % (self.symbol, n, 'spdf'[l])
            else:
                id = '%s-%s%d' % (self.symbol, 'spdf'[l], -n)
            self.id_j.append(id)
        setup.id_j = self.id_j

        if write_xml:
            setup.write_xml()
        return setup

    def diagonalize(self, h):
        ng = 350
        t = self.text
        t()
        t('Diagonalizing with gridspacing h=%.3f' % h)
        R = h * np.arange(1, ng + 1)
        G = (self.N * R / (self.beta + R) + 0.5).astype(int)
        G = np.clip(G, 1, self.N - 2)
        R1 = np.take(self.r, G - 1)
        R2 = np.take(self.r, G)
        R3 = np.take(self.r, G + 1)
        x1 = (R - R2) * (R - R3) / (R1 - R2) / (R1 - R3)
        x2 = (R - R1) * (R - R3) / (R2 - R1) / (R2 - R3)
        x3 = (R - R1) * (R - R2) / (R3 - R1) / (R3 - R2)
        def interpolate(f):
            f1 = np.take(f, G - 1)
            f2 = np.take(f, G)
            f3 = np.take(f, G + 1)
            return f1 * x1 + f2 * x2 + f3 * x3
        vt = interpolate(self.vt)
        t()
        t('state   all-electron     PAW')
        t('-------------------------------')
        for l in range(4):
            if l <= self.lmax:
                q_n = np.array([interpolate(q) for q in self.q_ln[l]])
                H = np.dot(np.transpose(q_n),
                           np.dot(self.dH_lnn[l], q_n)) * h
                S = np.dot(np.transpose(q_n),
                           np.dot(self.dO_lnn[l], q_n)) * h
            else:
                H = np.zeros((ng, ng))
                S = np.zeros((ng, ng))
            H.ravel()[::ng + 1] += vt + 1.0 / h**2 + l * (l + 1) / 2.0 / R**2
            H.ravel()[1::ng + 1] -= 0.5 / h**2
            H.ravel()[ng::ng + 1] -= 0.5 / h**2
            S.ravel()[::ng + 1] += 1.0
            e_n = np.zeros(ng)
            general_diagonalize(H, e_n, S)
            ePAW = e_n[0]
            if l <= self.lmax and self.n_ln[l][0] > 0:
                eAE = self.e_ln[l][0]
                t('%d%s:   %12.6f %12.6f' % (self.n_ln[l][0],
                                             'spdf'[l], eAE, ePAW), end='')
                if abs(eAE - ePAW) > 0.014:
                    t('  GHOST-STATE!')
                    self.ghost = True
                else:
                    t()
            else:
                t('*%s:                %12.6f' % ('spdf'[l], ePAW), end='')
                if ePAW < self.emax:
                    t('  GHOST-STATE!')
                    self.ghost = True
                else:
                    t()
        t('-------------------------------')

    def integrate(self, l, vt, e, gld, q=None):
        r = self.r[1:]
        dr = self.dr[1:]
        s = np.zeros(self.N)

        c0 = 0.5 * l * (l + 1) / r**2
        c1 = -0.5 * self.d2gdr2[1:]
        c2 = -0.5 * dr**-2

        fp = c2 + 0.5 * c1
        fm = c2 - 0.5 * c1
        f0 = c0 - 2 * c2

        f0 += vt[1:] - e
        if q is None:
            s[1] = r[1]**(l + 1)
            for g in range(gld):
                s[g + 2] = (-fm[g] * s[g] - f0[g] * s[g + 1]) / fp[g]
            return s

        s[1] = q[1] / (vt[0] - e)
        for g in range(gld):
            s[g + 2] = (q[g + 1] - fm[g] * s[g] - f0[g] * s[g + 1]) / fp[g]
        return s

    def write_xml(self, vl_j, vn_j, vf_j, ve_j, vu_j, vs_j, vq_j,
                  nc, nct, nt, Ekincore, X_p, ExxC, vbar,
                  tauc, tauct, extra_xc_data):
        raise DeprecationWarning('use gpaw/setup_data.py')
        xcname = self.xc.name
        if self.name is None:
            xml = open('%s.%s' % (self.symbol, xcname), 'w')
        else:
            xml = open('%s.%s.%s' % (self.symbol, self.name, xcname), 'w')
        if self.ghost:
            raise RuntimeError('Ghost!')

        print('<?xml version="1.0"?>', file=xml)
        print('<paw_setup version="0.6">', file=xml)

        name = atomic_names[self.Z].title()
        comment1 = name + ' setup for the Projector Augmented Wave method.'
        comment2 = 'Units: Hartree and Bohr radii.'
        comment2 += ' ' * (len(comment1) - len(comment2))
        print('  <!--', comment1, '-->', file=xml)
        print('  <!--', comment2, '-->', file=xml)

        print(('  <atom symbol="%s" Z="%d" core="%.1f" valence="%d"/>'
                       % (self.symbol, self.Z, self.Nc, self.Nv)), file=xml)
        if self.xcname == 'LDA':
            type = 'LDA'
            name = 'PW'
        else:
            type = 'GGA'
            name = self.xcname

        print('  <xc_functional type="%s" name="%s"/>' % (type, name), file=xml)
        if self.scalarrel:
            type = 'scalar-relativistic'
        else:
            type = 'non-relativistic'

        print('  <generator type="%s" name="gpaw-%s">' % \
              (type, version), file=xml)
        print('    Frozen core:', self.core or 'none', file=xml)
        print('  </generator>', file=xml)

        print('  <ae_energy kinetic="%f" xc="%f"' % \
              (self.Ekin, self.Exc), file=xml)
        print('             electrostatic="%f" total="%f"/>' % \
              (self.Epot, self.Ekin + self.Exc + self.Epot), file=xml)

        print('  <core_energy kinetic="%f"/>' % Ekincore, file=xml)

        print('  <valence_states>', file=xml)
        ids = []
        line1 = '    <state n="%d" l="%d" f=%s rc="%5.3f" e="%8.5f" id="%s"/>'
        line2 = '    <state       l="%d"        rc="%5.3f" e="%8.5f" id="%s"/>'
        for l, n, f, e in zip(vl_j, vn_j, vf_j, ve_j):
            if n > 0:
                f = '%-4s' % ('"%d"' % f)
                id = '%s-%d%s' % (self.symbol, n, 'spdf'[l])
                print(line1 % (n, l, f, self.rcut_l[l], e, id), file=xml)
            else:
                id = '%s-%s%d' % (self.symbol, 'spdf'[l], -n)
                print(line2 % (l, self.rcut_l[l], e, id), file=xml)
            ids.append(id)
        print('  </valence_states>', file=xml)

        print(('  <radial_grid eq="r=a*i/(n-i)" a="%f" n="%d" ' +
                       'istart="0" iend="%d" id="g1"/>') % \
                       (self.beta, self.N, self.N - 1), file=xml)

        rcgauss = self.rcutcomp / sqrt(self.gamma)
        print(('  <shape_function type="gauss" rc="%.12e"/>' %
                       rcgauss), file=xml)

        r = self.r

        if self.jcorehole is not None:
            print("self.jcorehole", self.jcorehole)
            print((('  <core_hole_state state="%d%s" ' +
                           'removed="%.1f" eig="%.8f" ekin="%.8f">') %
                           (self.ncorehole, 'spdf'[self.lcorehole],
                            self.fcorehole,
                            self.e_j[self.jcorehole],self.Ekincorehole)), file=xml)
            #print 'normalized?', np.dot(self.dr, self.u_j[self.jcorehole]**2)
            p = self.u_j[self.jcorehole].copy()
            p[1:] /= r[1:]
            if self.l_j[self.jcorehole] == 0:
                p[0] = (p[2] +
                        (p[1] - p[2]) * (r[0] - r[2]) / (r[1] - r[2]))
            for x in p:
                print('%16.12e' % x, end=' ', file=xml)
            print('\n  </core_hole_state>', file=xml)

        for name, a in [('ae_core_density', nc),
                        ('pseudo_core_density', nct),
                        ('pseudo_valence_density', nt - nct),
                        ('zero_potential', vbar),
                        ('ae_core_kinetic_energy_density',tauc),
                        ('pseudo_core_kinetic_energy_density',tauct)]:
            print('  <%s grid="g1">\n    ' % name, end=' ', file=xml)
            for x in a * sqrt(4 * pi):
                print('%16.12e' % x, end=' ', file=xml)
            print('\n  </%s>' % name, file=xml)

        # Print xc-specific data to setup file (used so for KLI and GLLB)
        for name, a in extra_xc_data.iteritems():
            newname = 'GLLB_'+name
            print('  <%s grid="g1">\n    ' % newname, end=' ', file=xml)
            for x in a:
                print('%16.12e' % x, end=' ', file=xml)
            print('\n  </%s>' % newname, file=xml)

        for l, u, s, q, in zip(vl_j, vu_j, vs_j, vq_j):
            id = ids.pop(0)
            for name, a in [('ae_partial_wave', u),
                            ('pseudo_partial_wave', s),
                            ('projector_function', q)]:
                print(('  <%s state="%s" grid="g1">\n    ' %
                               (name, id)), end=' ', file=xml)
                p = a.copy()
                p[1:] /= r[1:]
                if l == 0:
                    # XXXXX go to higher order!!!!!
                    p[0] = (p[2] +
                            (p[1] - p[2]) * (r[0] - r[2]) / (r[1] - r[2]))
                for x in p:
                    print('%16.12e' % x, end=' ', file=xml)
                print('\n  </%s>' % name, file=xml)

        print('  <kinetic_energy_differences>', end=' ', file=xml)
        nj = len(self.dK_jj)
        for j1 in range(nj):
            print('\n    ', end=' ', file=xml)
            for j2 in range(nj):
                print('%16.12e' % self.dK_jj[j1, j2], end=' ', file=xml)
        print('\n  </kinetic_energy_differences>', file=xml)

        if X_p is not None:
            print('  <exact_exchange_X_matrix>\n    ', end=' ', file=xml)
            for x in X_p:
                print('%16.12e' % x, end=' ', file=xml)
            print('\n  </exact_exchange_X_matrix>', file=xml)

            print('  <exact_exchange core-core="%f"/>' % ExxC, file=xml)

        print('</paw_setup>', file=xml)

def construct_smooth_wavefunction(u, l, gc, r, s):
    # Do a linear regression to a wave function
    # s = a + br^2 + cr^4 + dr^6, such that
    # the fitting is as good as possible in region gc-2:gc+2
    A = np.ones((4, 4))
    A[:, 0] = 1.0
    A[:, 1] = r[gc - 2:gc + 2]**2
    A[:, 2] = A[:, 1]**2
    A[:, 3] = A[:, 1] * A[:, 2]
    a = u[gc - 2:gc + 2] / r[gc - 2:gc + 2]**(l + 1)
    a = solve(A, a)
    r1 = r[:gc]
    r2 = r1**2
    rl1 = r1**(l + 1)
    y = a[0] + r2 * (a[1] + r2 * (a[2] + r2 * (a[3])))
    s[:gc] = rl1 * y


if __name__ == '__main__':
    import os
    from gpaw.atom.basis import BasisMaker
    from gpaw.atom.configurations import parameters

    for xcname in ['LDA', 'PBE', 'RPBE', 'revPBE', 'GLLBSC']:
        for symbol, par in parameters.items():
            filename = symbol + '.' + xcname
            if os.path.isfile(filename) or os.path.isfile(filename + '.gz'):
                continue
            g = Generator(symbol, xcname, scalarrel=True, nofiles=True)
            g.run(exx=True, logderiv=False, use_restart_file=False, **par)

            if xcname == 'PBE':
                bm = BasisMaker(g, name='dzp', run=False)
                basis = bm.generate()
                basis.write_xml()

