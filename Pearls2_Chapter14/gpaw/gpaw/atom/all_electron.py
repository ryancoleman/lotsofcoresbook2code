# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""
Atomic Density Functional Theory
"""

from __future__ import print_function
from math import pi, sqrt, log
import tempfile
import pickle
import sys
import os

import numpy as np
from ase.data import atomic_names
from ase.utils import devnull

from gpaw.atom.configurations import configurations
from gpaw.atom.radialgd import AERadialGridDescriptor
from gpaw.xc import XC
from gpaw.utilities import hartree
from gpaw import ConvergenceError

tempdir = tempfile.gettempdir()


# fine-structure constant
alpha = 1 / 137.036


class AllElectron:
    """Object for doing an atomic DFT calculation."""

    def __init__(self, symbol, xcname='LDA', scalarrel=False,
                 corehole=None, configuration=None, nofiles=True,
                 txt='-', gpernode=150, orbital_free=False, tf_coeff=1.):
        """Do an atomic DFT calculation.

        Example::

          a = AllElectron('Fe')
          a.run()
        """
        
        if txt is None:
            txt = devnull
        elif txt == '-':
            txt = sys.stdout
        elif isinstance(txt, str):
            txt = open(txt, 'w')
        self.txt = txt

        self.symbol = symbol
        self.xcname = xcname
        self.scalarrel = scalarrel
        self.nofiles = nofiles

        # Get reference state:
        self.Z, nlfe_j = configurations[symbol]

        # Collect principal quantum numbers, angular momentum quantum
        # numbers, occupation numbers and eigenvalues (j is a combined
        # index for n and l):
        self.n_j = [n for n, l, f, e in nlfe_j]
        self.l_j = [l for n, l, f, e in nlfe_j]
        self.f_j = [f for n, l, f, e in nlfe_j]
        self.e_j = [e for n, l, f, e in nlfe_j]

        if configuration is not None:
            j = 0
            for conf in configuration.split(','):
                if conf[0].isdigit():
                    n = int(conf[0])
                    l = 'spdf'.find(conf[1])
                    if len(conf) == 2:
                        f = 1.0
                    else:
                        f = float(conf[2:])
                    #try:
                    assert n == self.n_j[j]
                    assert l == self.l_j[j]
                    self.f_j[j] = f
                    #except IndexError:
                    #    self.n_j.append(n)
                    #    self.l_j.append(l)
                    #    self.f_j.append(f)
                    #    self.e_j.append(self.e_j[-1])
                    j += 1
                else:
                    j += {'He': 1,
                          'Ne': 3,
                          'Ar': 5,
                          'Kr': 8,
                          'Xe': 11}[conf]

        maxnodes = max([n - l - 1 for n, l in zip(self.n_j, self.l_j)])
        self.N = (maxnodes + 1) * gpernode
        self.beta = 0.4

        self.orbital_free = orbital_free
        self.tf_coeff = tf_coeff

        if self.orbital_free:
            self.n_j = [1]
            self.l_j = [0]
            self.f_j = [self.Z]
            self.e_j = [self.e_j[-1]]
            
        t = self.text
        t()
        if scalarrel:
            t('Scalar-relativistic atomic ', end='')
        else:
            t('Atomic ', end='')
        t('%s calculation for %s (%s, Z=%d)' % (xcname, symbol,
                                                atomic_names[self.Z], self.Z))

        if corehole is not None:
            self.ncorehole, self.lcorehole, self.fcorehole = corehole

            # Find j for core hole and adjust occupation:
            for j in range(len(self.f_j)):
                if (self.n_j[j] == self.ncorehole and
                    self.l_j[j] == self.lcorehole):
                    assert self.f_j[j] == 2 * (2 * self.lcorehole + 1)
                    self.f_j[j] -= self.fcorehole
                    self.jcorehole = j
                    break

            coreholestate = '%d%s' % (self.ncorehole, 'spdf'[self.lcorehole])
            t('Core hole in %s state (%s occupation: %.1f)' % (
                coreholestate, coreholestate, self.f_j[self.jcorehole]))
        else:
            self.jcorehole = None
            self.fcorehole = 0

    def text(self, *args, **kwargs):
        self.txt.write(kwargs.get('sep', ' ').join([str(arg)
                                                    for arg in args]) +
                       kwargs.get('end', '\n'))

    def initialize_wave_functions(self):
        r = self.r
        dr = self.dr
        # Initialize with Slater function:
        for l, e, u in zip(self.l_j, self.e_j, self.u_j):
            if self.symbol in ['Hf', 'Ta', 'W', 'Re', 'Os',
                               'Ir', 'Pt', 'Au']:
                a = sqrt(-4.0 * e)
            else:
                a = sqrt(-2.0 * e)

            u[:] = r**(1 + l) * np.exp(-a * r)
            norm = np.dot(u**2, dr)
            u *= 1.0 / sqrt(norm)
            
    def run(self, use_restart_file=True):
        #     beta g
        # r = ------, g = 0, 1, ..., N - 1
        #     N - g
        #
        #        rN
        # g = --------
        #     beta + r

        t = self.text
        N = self.N
        beta = self.beta
        t(N, 'radial gridpoints.')
        self.rgd = AERadialGridDescriptor(beta / N, 1.0 / N, N)
        g = np.arange(N, dtype=float)
        self.r = self.rgd.r_g
        self.dr = self.rgd.dr_g
        self.d2gdr2 = self.rgd.d2gdr2()

        # Number of orbitals:
        nj = len(self.n_j)

        # Radial wave functions multiplied by radius:
        self.u_j = np.zeros((nj, self.N))

        # Effective potential multiplied by radius:
        self.vr = np.zeros(N)

        # Electron density:
        self.n = np.zeros(N)

        # Always spinpaired nspins=1
        self.xc = XC(self.xcname)

        # Initialize for non-local functionals
        if self.xc.type == 'GLLB':
            self.xc.pass_stuff_1d(self)
            self.xc.initialize_1d()
            
        n_j = self.n_j
        l_j = self.l_j
        f_j = self.f_j
        e_j = self.e_j
        
        Z = self.Z    # nuclear charge
        r = self.r    # radial coordinate
        dr = self.dr  # dr/dg
        n = self.n    # electron density
        vr = self.vr  # effective potential multiplied by r

        vHr = np.zeros(self.N)
        self.vXC = np.zeros(self.N)

        restartfile = '%s/%s.restart' % (tempdir, self.symbol)
        if self.xc.type == 'GLLB' or not use_restart_file:
            # Do not start from initial guess when doing
            # non local XC!
            # This is because we need wavefunctions as well
            # on the first iteration.
            fd = None
        else:
            try:
                fd = open(restartfile, 'r')
            except IOError:
                fd = None
            else:
                try:
                    n[:] = pickle.load(fd)
                except (ValueError, IndexError):
                    fd = None
                else:
                    norm = np.dot(n * r**2, dr) * 4 * pi
                    if abs(norm - sum(f_j)) > 0.01:
                        fd = None
                    else:
                        t('Using old density for initial guess.')
                        n *= sum(f_j) / norm

        if fd is None:
            self.initialize_wave_functions()
            n[:] = self.calculate_density()

        bar = '|------------------------------------------------|'
        t(bar)
        
        niter = 0
        nitermax = 117
        qOK = log(1e-10)
        mix = 0.4
        
        # orbital_free needs more iterations and coefficient
        if self.orbital_free:
            #qOK = log(1e-14)
            e_j[0] /= self.tf_coeff
            mix = 0.01
            nitermax = 1000
            
        vrold = None
        
        while True:
            # calculate hartree potential
            hartree(0, n * r * dr, r, vHr)

            # add potential from nuclear point charge (v = -Z / r)
            vHr -= Z

            # calculated exchange correlation potential and energy
            self.vXC[:] = 0.0

            if self.xc.type == 'GLLB':
                # Update the potential to self.vXC an the energy to self.Exc
                Exc = self.xc.get_xc_potential_and_energy_1d(self.vXC)
            else:
                Exc = self.xc.calculate_spherical(self.rgd,
                                                  n.reshape((1, -1)),
                                                  self.vXC.reshape((1, -1)))

            # calculate new total Kohn-Sham effective potential and
            # admix with old version
            vr[:] = (vHr + self.vXC * r) / self.tf_coeff

            if niter > 0:
                vr[:] = mix * vr + (1 - mix) * vrold
            vrold = vr.copy()

            # solve Kohn-Sham equation and determine the density change
            self.solve()
            dn = self.calculate_density() - n
            n += dn

            # estimate error from the square of the density change integrated
            q = log(np.sum((r * dn)**2))

            # print progress bar
            if niter == 0:
                q0 = q
                b0 = 0
            else:
                b = int((q0 - q) / (q0 - qOK) * 50)
                if b > b0:
                    self.txt.write(bar[b0:min(b, 50)])
                    self.txt.flush()
                    b0 = b

            # check if converged and break loop if so
            if q < qOK:
                self.txt.write(bar[b0:])
                self.txt.flush()
                break

            niter += 1
            if niter > nitermax:
                raise RuntimeError('Did not converge!')

        tau = self.calculate_kinetic_energy_density()

        t()
        t('Converged in %d iteration%s.' % (niter, 's'[:niter != 1]))

        try:
            fd = open(restartfile, 'w')
        except IOError:
            pass
        else:
            pickle.dump(n, fd)
            try:
                os.chmod(restartfile, 0666)
            except OSError:
                pass

        Ekin = 0
        if self.orbital_free:
            e_j[0] *= self.tf_coeff
            vr *= self.tf_coeff
        
        for f, e in zip(f_j, e_j):
            Ekin += f * e

        Epot = 2 * pi * np.dot(n * r * (vHr - Z), dr)
        Ekin += -4 * pi * np.dot(n * vr * r, dr)

        t()
        t('Energy contributions:')
        t('-------------------------')
        t('Kinetic:   %+13.6f' % Ekin)
        t('XC:        %+13.6f' % Exc)
        t('Potential: %+13.6f' % Epot)
        t('-------------------------')
        t('Total:     %+13.6f' % (Ekin + Exc + Epot))
        self.ETotal = Ekin + Exc + Epot
        t()

        t('state      eigenvalue         ekin         rmax')
        t('-----------------------------------------------')
        for m, l, f, e, u in zip(n_j, l_j, f_j, e_j, self.u_j):
            # Find kinetic energy:
            k = e - np.sum((np.where(abs(u) < 1e-160, 0, u)**2 *  # XXXNumeric!
                            vr * dr)[1:] / r[1:])

            # Find outermost maximum:
            g = self.N - 4
            while u[g - 1] >= u[g]:
                g -= 1
            x = r[g - 1:g + 2]
            y = u[g - 1:g + 2]
            A = np.transpose(np.array([x**i for i in range(3)]))
            c, b, a = np.linalg.solve(A, y)
            assert a < 0.0
            rmax = -0.5 * b / a

            s = 'spdf'[l]
            t('%d%s^%-4.1f: %12.6f %12.6f %12.3f' % (m, s, f, e, k, rmax))
        t('-----------------------------------------------')
        t('(units: Bohr and Hartree)')

        for m, l, u in zip(n_j, l_j, self.u_j):
            self.write(u, 'ae', n=m, l=l)

        self.write(n, 'n')
        self.write(vr, 'vr')
        self.write(vHr, 'vHr')
        self.write(self.vXC, 'vXC')
        self.write(tau, 'tau')

        self.Ekin = Ekin
        self.Epot = Epot
        self.Exc = Exc

    def write(self, array, name=None, n=None, l=None):
        if self.nofiles:
            return

        if name:
            name = self.symbol + '.' + name
        else:
            name = self.symbol

        if l is not None:
            assert n is not None
            if n > 0:
                # Bound state:
                name += '.%d%s' % (n, 'spdf'[l])
            else:
                name += '.x%d%s' % (-n, 'spdf'[l])

        f = open(name, 'w')
        for r, a in zip(self.r, array):
            print(r, a, file=f)

    def calculate_density(self):
        """Return the electron charge density divided by 4 pi"""
        n = np.dot(self.f_j,
                   np.where(abs(self.u_j) < 1e-160, 0,
                            self.u_j)**2) / (4 * pi)
        n[1:] /= self.r[1:]**2
        n[0] = n[1]
        return n

    def calculate_kinetic_energy_density(self):
        """Return the kinetic energy density"""
        return self.radial_kinetic_energy_density(self.f_j, self.l_j, self.u_j)

    def radial_kinetic_energy_density(self, f_j, l_j, u_j):
        """Kinetic energy density from a restricted set of wf's
        """
        shape = np.shape(u_j[0])
        dudr = np.zeros(shape)
        tau = np.zeros(shape)
        for f, l, u in zip(f_j, l_j, u_j):
            self.rgd.derivative(u, dudr)
            # contribution from angular derivatives
            if l > 0:
                tau += f * l * (l + 1) * np.where(abs(u) < 1e-160, 0, u)**2
            # contribution from radial derivatives
            dudr = u - self.r * dudr
            tau += f * np.where(abs(dudr) < 1e-160, 0, dudr)**2
        tau[1:] /= self.r[1:]**4
        tau[0] = tau[1]

        return 0.5 * tau / (4 * pi)

    def calculate_kinetic_energy_density2(self):
        """Return the kinetic energy density
        calculation over R(r)=u(r)/r
        slower convergence with # of radial grid points for
        Ekin of H than radial_kinetic_energy_density
        """

        shape = self.u_j.shape[1]
        R = np.zeros(shape)
        dRdr = np.zeros(shape)
        tau = np.zeros(shape)
        for f, l, u in zip(self.f_j, self.l_j, self.u_j):
            R[1:] = u[1:] / self.r[1:]
            if l == 0:
                # estimate value at origin by Taylor series to first order
                d1 = self.r[1]
                d2 = self.r[2]
                R[0] = .5 * (R[1] + R[2] + (R[1] - R[2]) *
                             (d1 + d2) / (d2 - d1))
            else:
                R[0] = 0
            self.rgd.derivative(R, dRdr)
            # contribution from radial derivatives
            tau += f * np.where(abs(dRdr) < 1e-160, 0, dRdr)**2
            # contribution from angular derivatives
            if l > 0:
                R[1:] = R[1:] / self.r[1:]
                if l == 1:
                    R[0] = R[1]
                else:
                    R[0] = 0
                tau += f * l * (l + 1) * np.where(abs(R) < 1e-160, 0, R)**2

        return 0.5 * tau / (4 * pi)

    def solve(self):
        """Solve the Schrodinger equation

        ::

             2
            d u     1  dv  du   u     l(l + 1)
          - --- - ---- -- (-- - -) + [-------- + 2M(v - e)] u = 0,
              2      2 dr  dr   r         2
            dr    2Mc                    r


        where the relativistic mass::

                   1
          M = 1 - --- (v - e)
                    2
                  2c

        and the fine-structure constant alpha = 1/c = 1/137.036
        is set to zero for non-scalar-relativistic calculations.

        On the logaritmic radial grids defined by::

              beta g
          r = ------,  g = 0, 1, ..., N - 1
              N - g

                 rN
          g = --------, r = [0; oo[
              beta + r

        the Schrodinger equation becomes::

           2
          d u      du
          --- c  + -- c  + u c  = 0
            2  2   dg  1      0
          dg

        with the vectors c0, c1, and c2  defined by::

                 2 dg 2
          c  = -r (--)
           2       dr

                  2         2
                 d g  2    r   dg dv
          c  = - --- r  - ---- -- --
           1       2         2 dr dr
                 dr       2Mc

                                    2    r   dv
          c  = l(l + 1) + 2M(v - e)r  + ---- --
           0                               2 dr
                                        2Mc
        """
        r = self.r
        dr = self.dr
        vr = self.vr

        c2 = -(r / dr)**2
        c10 = -self.d2gdr2 * r**2  # first part of c1 vector

        if self.scalarrel:
            self.r2dvdr = np.zeros(self.N)
            self.rgd.derivative(vr, self.r2dvdr)
            self.r2dvdr *= r
            self.r2dvdr -= vr
        else:
            self.r2dvdr = None

        # solve for each quantum state separately
        for j, (n, l, e, u) in enumerate(zip(self.n_j, self.l_j,
                                             self.e_j, self.u_j)):
            nodes = n - l - 1  # analytically expected number of nodes
            delta = -0.2 * e
            nn, A = shoot(u, l, vr, e, self.r2dvdr, r, dr, c10, c2,
                          self.scalarrel)
            # adjust eigenenergy until u has the correct number of nodes
            while nn != nodes:
                diff = cmp(nn, nodes)
                while diff == cmp(nn, nodes):
                    e -= diff * delta
                    nn, A = shoot(u, l, vr, e, self.r2dvdr, r, dr, c10, c2,
                                  self.scalarrel)
                delta /= 2

            # adjust eigenenergy until u is smooth at the turning point
            de = 1.0
            while abs(de) > 1e-9:
                norm = np.dot(np.where(abs(u) < 1e-160, 0, u)**2, dr)
                u *= 1.0 / sqrt(norm)
                de = 0.5 * A / norm
                x = abs(de / e)
                if x > 0.1:
                    de *= 0.1 / x
                e -= de
                assert e < 0.0
                nn, A = shoot(u, l, vr, e, self.r2dvdr, r, dr, c10, c2,
                              self.scalarrel)
            self.e_j[j] = e
            u *= 1.0 / sqrt(np.dot(np.where(abs(u) < 1e-160, 0, u)**2, dr))

    def solve_confined(self, j, rc, vconf=None):
        """Solve the Schroedinger equation in a confinement potential.
        
        Solves the Schroedinger equation like the solve method, but with a
        number of differences.  Before invoking this method, run solve() to
        get initial guesses.

        Parameters:
            j: solves only for the state given by j
            rc: solution cutoff. Solution will be zero outside this.
            vconf: added to the potential (use this as confinement potential)

        Returns: a tuple containing the solution u and its energy e.

        Unlike the solve method, this method will not alter any attributes of
        this object.
        """
        r = self.r
        dr = self.dr
        vr = self.vr.copy()
        if vconf is not None:
            vr += vconf * r

        c2 = -(r / dr)**2
        c10 = -self.d2gdr2 * r**2  # first part of c1 vector

        if j is None:
            n, l, e, u = 3, 2, -0.15, self.u_j[-1].copy()
        else:
            n = self.n_j[j]
            l = self.l_j[j]
            e = self.e_j[j]
            u = self.u_j[j].copy()
            
        nn, A = shoot_confined(u, l, vr, e, self.r2dvdr, r, dr, c10, c2,
                               self.scalarrel, rc=rc, beta=self.beta)
        assert nn == n - l - 1  # run() should have been called already
        
        # adjust eigenenergy until u is smooth at the turning point
        de = 1.0
        while abs(de) > 1e-9:
            norm = np.dot(np.where(abs(u) < 1e-160, 0, u)**2, dr)
            u *= 1.0 / sqrt(norm)
            de = 0.5 * A / norm
            x = abs(de / e)
            if x > 0.1:
                de *= 0.1 / x
            e -= de
            assert e < 0.0

            nn, A = shoot_confined(u, l, vr, e, self.r2dvdr, r, dr, c10, c2,
                                   self.scalarrel, rc=rc, beta=self.beta)
        u *= 1.0 / sqrt(np.dot(np.where(abs(u) < 1e-160, 0, u)**2, dr))
        return u, e

    def kin(self, l, u, e=None):  # XXX move to Generator
        r = self.r[1:]
        dr = self.dr[1:]

        c0 = 0.5 * l * (l + 1) / r**2
        c1 = -0.5 * self.d2gdr2[1:]
        c2 = -0.5 * dr**-2

        if e is not None and self.scalarrel:
            x = 0.5 * alpha**2
            Mr = r * (1.0 + x * e) - x * self.vr[1:]
            c0 += ((Mr - r) * (self.vr[1:] - e * r) +
                   0.5 * x * self.r2dvdr[1:] / Mr) / r**2
            c1 -= 0.5 * x * self.r2dvdr[1:] / (Mr * dr * r)

        fp = c2 + 0.5 * c1
        fm = c2 - 0.5 * c1
        f0 = c0 - 2 * c2
        kr = np.zeros(self.N)
        kr[1:] = f0 * u[1:] + fm * u[:-1]
        kr[1:-1] += fp[:-1] * u[2:]
        kr[0] = 0.0
        return kr

    def r2g(self, r):
        """Convert radius to index of the radial grid."""
        return int(r * self.N / (self.beta + r))

    def get_confinement_potential(self, alpha, ri, rc):
        """Create a smooth confinement potential.
        
        Returns a (potential) function which is zero inside the radius ri
        and goes to infinity smoothly at rc, after which point it is nan.
        The potential is given by::

                   alpha         /   rc - ri \
          V(r) = --------   exp ( - --------- )   for   ri < r < rc
                  rc - r         \    r - ri /

        """
        i_ri = self.r2g(ri)
        i_rc = self.r2g(rc)
        if self.r[i_rc] == rc:
            # Avoid division by zero in the odd case that rc coincides
            # exactly with a grid point (which actually happens sometimes)
            i_rc -= 1

        potential = np.zeros(np.shape(self.r))
        r = self.r[i_ri + 1:i_rc + 1]
        exponent = - (rc - ri) / (r - ri)
        denom = rc - r
        value = np.exp(exponent) / denom
        potential[i_ri + 1:i_rc + 1] = value
        potential[i_rc + 1:] = np.inf

        return alpha * potential


def shoot(u, l, vr, e, r2dvdr, r, dr, c10, c2, scalarrel=False, gmax=None):
    """n, A = shoot(u, l, vr, e, ...)

    For guessed trial eigenenergy e, integrate the radial Schrodinger
    equation::

          2
         d u      du
         --- c  + -- c  + u c  = 0
           2  2   dg  1      0
         dg

               2 dg 2
        c  = -r (--)
         2       dr

                2         2
               d g  2    r   dg dv
        c  = - --- r  - ---- -- --
         1       2         2 dr dr
               dr       2Mc

                                  2    r   dv
        c  = l(l + 1) + 2M(v - e)r  + ---- --
         0                               2 dr
                                      2Mc

    The resulting wavefunction is returned in input vector u.
    The number of nodes of u is returned in attribute n.
    Returned attribute A, is a measure of the size of the derivative
    discontinuity at the classical turning point.
    The trial energy e is correct if A is zero and n is the correct number
    of nodes."""

    if scalarrel:
        x = 0.5 * alpha**2  # x = 1 / (2c^2)
        Mr = r * (1.0 + x * e) - x * vr
    else:
        Mr = r
    c0 = l * (l + 1) + 2 * Mr * (vr - e * r)
    if gmax is None and np.alltrue(c0 > 0):
        raise ConvergenceError('Bad initial electron density guess!')

    c1 = c10
    if scalarrel:
        c0 += x * r2dvdr / Mr
        c1 = c10 - x * r * r2dvdr / (Mr * dr)

    # vectors needed for numeric integration of diff. equation
    fm = 0.5 * c1 - c2
    fp = 0.5 * c1 + c2
    f0 = c0 - 2 * c2

    if gmax is None:
        # set boundary conditions at r -> oo (u(oo) = 0 is implicit)
        u[-1] = 1.0

        # perform backwards integration from infinity to the turning point
        g = len(u) - 2
        u[-2] = u[-1] * f0[-1] / fm[-1]
        while c0[g] > 0.0:  # this defines the classical turning point
            u[g - 1] = (f0[g] * u[g] + fp[g] * u[g + 1]) / fm[g]
            if u[g - 1] < 0.0:
                # There should't be a node here!  Use a more negative
                # eigenvalue:
                print('!!!!!!', end=' ')
                return 100, None
            if u[g - 1] > 1e100:
                u *= 1e-100
            g -= 1

        # stored values of the wavefunction and the first derivative
        # at the turning point
        gtp = g + 1
        utp = u[gtp]
        dudrplus = 0.5 * (u[gtp + 1] - u[gtp - 1]) / dr[gtp]
    else:
        gtp = gmax

    # set boundary conditions at r -> 0
    u[0] = 0.0
    u[1] = 1.0

    # perform forward integration from zero to the turning point
    g = 1
    nodes = 0
    while g <= gtp:  # integrate one step further than gtp
                     # (such that dudr is defined in gtp)
        u[g + 1] = (fm[g] * u[g - 1] - f0[g] * u[g]) / fp[g]
        if u[g + 1] * u[g] < 0:
            nodes += 1
        g += 1
    if gmax is not None:
        return

    # scale first part of wavefunction, such that it is continuous at gtp
    u[:gtp + 2] *= utp / u[gtp]

    # determine size of the derivative discontinuity at gtp
    dudrminus = 0.5 * (u[gtp + 1] - u[gtp - 1]) / dr[gtp]
    A = (dudrplus - dudrminus) * utp

    return nodes, A


def shoot_confined(u, l, vr, e, r2dvdr, r, dr, c10, c2, scalarrel=False,
                   gmax=None, rc=10., beta=7.):
    """This method is used by the solve_confined method."""
    # XXX much of this is pasted from the ordinary shoot method

    if scalarrel:
        x = 0.5 * alpha**2  # x = 1 / (2c^2)
        Mr = r * (1.0 + x * e) - x * vr
    else:
        Mr = r
    c0 = l * (l + 1) + 2 * Mr * (vr - e * r)
    if gmax is None and np.alltrue(c0 > 0):
        print("""
Problem with initial electron density guess!  Try to run the program
with the '-n' option (non-scalar-relativistic calculation) and then
try again without the '-n' option (this will generate a good initial
guess for the density).
""")
        raise SystemExit
    c1 = c10
    if scalarrel:
        c0 += x * r2dvdr / Mr
        c1 = c10 - x * r * r2dvdr / (Mr * dr)

    # vectors needed for numeric integration of diff. equation
    fm = 0.5 * c1 - c2
    fp = 0.5 * c1 + c2
    f0 = c0 - 2 * c2

    if gmax is None:
        gcut = int(rc * len(r) / (beta + rc))
        # set boundary conditions at r -> oo (u(oo) = 0 is implicit)
        u[gcut - 1] = 1.
        u[gcut:] = 0.

        # perform backwards integration from infinity to the turning point
        g = gcut - 2
        u[g] = u[g + 1] * f0[g + 1] / fm[g + 1]
        
        while c0[g] > 0.0:  # this defines the classical turning point
            u[g - 1] = (f0[g] * u[g] + fp[g] * u[g + 1]) / fm[g]
            if u[g - 1] < 0.0:
                # There should't be a node here!  Use a more negative
                # eigenvalue:
                print('!!!!!!', end=' ')
                return 100, None
            if u[g - 1] > 1e100:
                u *= 1e-100
            g -= 1

        # stored values of the wavefunction and the first derivative
        # at the turning point
        gtp = g + 1
        utp = u[gtp]
        dudrplus = 0.5 * (u[gtp + 1] - u[gtp - 1]) / dr[gtp]
    else:
        gtp = gmax

    # set boundary conditions at r -> 0
    u[0] = 0.0
    u[1] = 1.0

    # perform forward integration from zero to the turning point
    g = 1
    nodes = 0
    while g <= gtp:  # integrate one step further than gtp
                     # (such that dudr is defined in gtp)
        u[g + 1] = (fm[g] * u[g - 1] - f0[g] * u[g]) / fp[g]
        if u[g + 1] * u[g] < 0:
            nodes += 1
        g += 1
    if gmax is not None:
        return

    # scale first part of wavefunction, such that it is continuous at gtp
    u[:gtp + 2] *= utp / u[gtp]

    # determine size of the derivative discontinuity at gtp
    dudrminus = 0.5 * (u[gtp + 1] - u[gtp - 1]) / dr[gtp]
    A = (dudrplus - dudrminus) * utp

    return nodes, A

if __name__ == '__main__':
    a = AllElectron('Cu', scalarrel=True)
    a.run()
