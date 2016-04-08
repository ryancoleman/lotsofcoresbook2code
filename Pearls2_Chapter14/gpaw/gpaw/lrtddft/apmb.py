"""Omega matrix for functionals with Hartree-Fock exchange.

"""
from __future__ import print_function
from math import sqrt

import numpy as np
from numpy.linalg import inv

from ase.units import Hartree
from ase.utils.timing import Timer

import gpaw.mpi as mpi
from gpaw.lrtddft.omega_matrix import OmegaMatrix
from gpaw.pair_density import PairDensity
from gpaw.utilities import pack
from gpaw.utilities.lapack import diagonalize, gemm, sqrt_matrix


class ApmB(OmegaMatrix):

    """Omega matrix for functionals with Hartree-Fock exchange.

    """

    def get_full(self):

        hybrid = ((self.xc is not None) and
                  hasattr(self.xc, 'hybrid') and
                  (self.xc.hybrid > 0.0))
        if self.fullkss.npspins < 2 and hybrid:
            raise RuntimeError('Does not work spin-unpolarized ' +
                               'with hybrids (use nspins=2)')

        self.paw.timer.start('ApmB RPA')
        self.ApB = self.Om
        self.AmB = self.get_rpa()
        self.paw.timer.stop()

        if self.xc is not None:
            self.paw.timer.start('ApmB XC')
            self.get_xc()  # inherited from OmegaMatrix
            self.paw.timer.stop()

    def get_rpa(self):
        """Calculate RPA and Hartree-fock part of the A+-B matrices."""

        # shorthands
        kss = self.fullkss
        finegrid = self.finegrid

        # calculate omega matrix
        nij = len(kss)
        print('RPAhyb', nij, 'transitions', file=self.txt)

        AmB = np.zeros((nij, nij))
        ApB = self.ApB

        # storage place for Coulomb integrals
        integrals = {}

        for ij in range(nij):
            print('RPAhyb kss[' + '%d' % ij + ']=', kss[ij], file=self.txt)

            timer = Timer()
            timer.start('init')
            timer2 = Timer()

            # smooth density including compensation charges
            timer2.start('with_compensation_charges 0')
            rhot_p = kss[ij].with_compensation_charges(
                finegrid is not 0)
            timer2.stop()

            # integrate with 1/|r_1-r_2|
            timer2.start('poisson')
            phit_p = np.zeros(rhot_p.shape, rhot_p.dtype)
            self.poisson.solve(phit_p, rhot_p, charge=None)
            timer2.stop()

            timer.stop()
            t0 = timer.get_time('init')
            timer.start(ij)

            if finegrid == 1:
                rhot = kss[ij].with_compensation_charges()
                phit = self.gd.zeros()
                self.restrict(phit_p, phit)
            else:
                phit = phit_p
                rhot = rhot_p

            for kq in range(ij, nij):
                if kq != ij:
                    # smooth density including compensation charges
                    timer2.start('kq with_compensation_charges')
                    rhot = kss[kq].with_compensation_charges(
                        finegrid is 2)
                    timer2.stop()
                pre = self.weight_Kijkq(ij, kq)

                timer2.start('integrate')
                I = self.Coulomb_integral_kss(kss[ij], kss[kq], phit, rhot)
                if kss[ij].spin == kss[kq].spin:
                    name = self.Coulomb_integral_name(kss[ij].i, kss[ij].j,
                                                      kss[kq].i, kss[kq].j,
                                                      kss[ij].spin)
                    integrals[name] = I
                ApB[ij, kq] = pre * I
                timer2.stop()

                if ij == kq:
                    epsij = kss[ij].get_energy() / kss[ij].get_weight()
                    AmB[ij, kq] += epsij
                    ApB[ij, kq] += epsij

            timer.stop()
# timer2.write()
            if ij < (nij - 1):
                print('RPAhyb estimated time left',
                      self.time_left(timer, t0, ij, nij), file=self.txt)

        # add HF parts and apply symmetry
        if hasattr(self.xc, 'hybrid'):
            weight = self.xc.hybrid
        else:
            weight = 0.0
        for ij in range(nij):
            print('HF kss[' + '%d' % ij + ']', file=self.txt)
            timer = Timer()
            timer.start('init')
            timer.stop()
            t0 = timer.get_time('init')
            timer.start(ij)

            i = kss[ij].i
            j = kss[ij].j
            s = kss[ij].spin
            for kq in range(ij, nij):
                if kss[ij].pspin == kss[kq].pspin:
                    k = kss[kq].i
                    q = kss[kq].j
                    ikjq = self.Coulomb_integral_ijkq(i, k, j, q, s, integrals)
                    iqkj = self.Coulomb_integral_ijkq(i, q, k, j, s, integrals)
                    ApB[ij, kq] -= weight * (ikjq + iqkj)
                    AmB[ij, kq] -= weight * (ikjq - iqkj)

                ApB[kq, ij] = ApB[ij, kq]
                AmB[kq, ij] = AmB[ij, kq]

            timer.stop()
            if ij < (nij - 1):
                print('HF estimated time left',
                      self.time_left(timer, t0, ij, nij), file=self.txt)

        return AmB

    def Coulomb_integral_name(self, i, j, k, l, spin):
        """return a unique name considering the Coulomb integral
        symmetry"""
        def ij_name(i, j):
            return str(max(i, j)) + ' ' + str(min(i, j))

        # maximal gives the first
        if max(i, j) >= max(k, l):
            base = ij_name(i, j) + ' ' + ij_name(k, l)
        else:
            base = ij_name(k, l) + ' ' + ij_name(i, j)
        return base + ' ' + str(spin)

    def Coulomb_integral_ijkq(self, i, j, k, q, spin, integrals):
        name = self.Coulomb_integral_name(i, j, k, q, spin)
        if name in integrals:
            return integrals[name]

        # create the Kohn-Sham singles
        kss_ij = PairDensity(self.paw)
        kss_ij.initialize(self.paw.wfs.kpt_u[spin], i, j)
        kss_kq = PairDensity(self.paw)
        kss_kq.initialize(self.paw.wfs.kpt_u[spin], k, q)

        rhot_p = kss_ij.with_compensation_charges(
            self.finegrid is not 0)
        phit_p = np.zeros(rhot_p.shape, rhot_p.dtype)
        self.poisson.solve(phit_p, rhot_p, charge=None)

        if self.finegrid == 1:
            phit = self.gd.zeros()
            self.restrict(phit_p, phit)
        else:
            phit = phit_p

        rhot = kss_kq.with_compensation_charges(
            self.finegrid is 2)

        integrals[name] = self.Coulomb_integral_kss(kss_ij, kss_kq,
                                                    phit, rhot)
        return integrals[name]

    def timestring(self, t):
        ti = int(t + .5)
        td = int(ti // 86400)
        st = ''
        if td > 0:
            st += '%d' % td + 'd'
            ti -= td * 86400
        th = int(ti // 3600)
        if th > 0:
            st += '%d' % th + 'h'
            ti -= th * 3600
        tm = int(ti // 60)
        if tm > 0:
            st += '%d' % tm + 'm'
            ti -= tm * 60
        st += '%d' % ti + 's'
        return st

    def map(self, istart=None, jend=None, energy_range=None):
        """Map A+B, A-B matrices according to constraints."""

        map, self.kss = self.get_map(istart, jend, energy_range)
        if map is None:
            ApB = self.ApB.copy()
            AmB = self.AmB.copy()
        else:
            nij = len(self.kss)
            ApB = np.empty((nij, nij))
            AmB = np.empty((nij, nij))
            for ij in range(nij):
                for kq in range(nij):
                    ApB[ij, kq] = self.ApB[map[ij], map[kq]]
                    AmB[ij, kq] = self.AmB[map[ij], map[kq]]

        return ApB, AmB

    def diagonalize(self, istart=None, jend=None, energy_range=None,
                    TDA=False):
        """Evaluate Eigenvectors and Eigenvalues"""

        ApB, AmB = self.map(istart, jend, energy_range)
        nij = len(self.kss)

        if TDA:
            # Tamm-Dancoff approximation (B=0)
            self.eigenvectors = 0.5 * (ApB + AmB)
            eigenvalues = np.zeros((nij))
            diagonalize(self.eigenvectors, eigenvalues)
            self.eigenvalues = eigenvalues ** 2
        else:
            # the occupation matrix
            C = np.empty((nij,))
            for ij in range(nij):
                C[ij] = 1. / self.kss[ij].fij

            S = C * inv(AmB) * C
            S = sqrt_matrix(inv(S).copy())

            # get Omega matrix
            M = np.zeros(ApB.shape)
            gemm(1.0, ApB, S, 0.0, M)
            self.eigenvectors = np.zeros(ApB.shape)
            gemm(1.0, S, M, 0.0, self.eigenvectors)

            self.eigenvalues = np.zeros((nij))
            diagonalize(self.eigenvectors, self.eigenvalues)

    def read(self, filename=None, fh=None):
        """Read myself from a file"""
        if mpi.rank == mpi.MASTER:
            if fh is None:
                f = open(filename, 'r')
            else:
                f = fh

            f.readline()
            nij = int(f.readline())
            ApB = np.zeros((nij, nij))
            for ij in range(nij):
                l = f.readline().split()
                for kq in range(ij, nij):
                    ApB[ij, kq] = float(l[kq - ij])
                    ApB[kq, ij] = ApB[ij, kq]
            self.ApB = ApB

            f.readline()
            nij = int(f.readline())
            AmB = np.zeros((nij, nij))
            for ij in range(nij):
                l = f.readline().split()
                for kq in range(ij, nij):
                    AmB[ij, kq] = float(l[kq - ij])
                    AmB[kq, ij] = AmB[ij, kq]
            self.AmB = AmB

            if fh is None:
                f.close()

    def weight_Kijkq(self, ij, kq):
        """weight for the coupling matrix terms"""
        return 2.

    def write(self, filename=None, fh=None):
        """Write current state to a file."""
        if mpi.rank == mpi.MASTER:
            if fh is None:
                f = open(filename, 'w')
            else:
                f = fh

            f.write('# A+B\n')
            nij = len(self.fullkss)
            f.write('%d\n' % nij)
            for ij in range(nij):
                for kq in range(ij, nij):
                    f.write(' %g' % self.ApB[ij, kq])
                f.write('\n')

            f.write('# A-B\n')
            nij = len(self.fullkss)
            f.write('%d\n' % nij)
            for ij in range(nij):
                for kq in range(ij, nij):
                    f.write(' %g' % self.AmB[ij, kq])
                f.write('\n')

            if fh is None:
                f.close()

    def __str__(self):
        string = '<ApmB> '
        if hasattr(self, 'eigenvalues'):
            string += 'dimension ' + ('%d' % len(self.eigenvalues))
            string += '\neigenvalues: '
            for ev in self.eigenvalues:
                string += ' ' + ('%f' % (sqrt(ev) * Hartree))
        return string
