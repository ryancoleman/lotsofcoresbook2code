from __future__ import print_function
import sys
from math import sqrt
import numpy as np
import gpaw.mpi as mpi
MASTER = mpi.MASTER

from ase.parallel import paropen
from ase.units import Hartree
from ase.utils.timing import Timer

from gpaw import debug
import gpaw.mpi as mpi
from gpaw.poisson import PoissonSolver
from gpaw.output import get_txt
from gpaw.lrtddft.excitation import Excitation, ExcitationList
from gpaw.lrtddft.kssingle import KSSingles
from gpaw.transformers import Transformer
from gpaw.utilities import pack, pack2, packed_index
from gpaw.utilities.lapack import diagonalize
from gpaw.xc import XC

import time

"""This module defines a Omega Matrix class."""


class OmegaMatrix:

    """
    Omega matrix in Casidas linear response formalism

    Parameters
      - calculator: the calculator object the ground state calculation
      - kss: the Kohn-Sham singles object
      - xc: the exchange correlation approx. to use
      - derivativeLevel: which level i of d^i Exc/dn^i to use
      - numscale: numeric epsilon for derivativeLevel=0,1
      - filehandle: the oject can be read from a filehandle
      - txt: output stream or file name
      - finegrid: level of fine grid to use. 0: nothing, 1 for poisson only,
        2 everything on the fine grid
    """

    def __init__(self,
                 calculator=None,
                 kss=None,
                 xc=None,
                 derivativeLevel=None,
                 numscale=0.001,
                 filehandle=None,
                 txt=None,
                 finegrid=2,
                 eh_comm=None,
                 ):

        if not txt and calculator:
            txt = calculator.txt
        self.txt = get_txt(txt, mpi.rank)

        if eh_comm is None:
            eh_comm = mpi.serial_comm

        self.eh_comm = eh_comm

        if filehandle is not None:
            self.kss = kss
            self.read(fh=filehandle)
            return None

        self.fullkss = kss
        self.finegrid = finegrid

        if calculator is None:
            return

        self.paw = calculator
        wfs = self.paw.wfs

        # handle different grid possibilities
        self.restrict = None
        # self.poisson = PoissonSolver(nn=self.paw.hamiltonian.poisson.nn)
        self.poisson = calculator.hamiltonian.poisson
        if finegrid:
            self.poisson.set_grid_descriptor(self.paw.density.finegd)
            self.poisson.initialize()

            self.gd = self.paw.density.finegd
            if finegrid == 1:
                self.gd = wfs.gd
        else:
            self.poisson.set_grid_descriptor(wfs.gd)
            self.poisson.initialize()
            self.gd = wfs.gd
        self.restrict = Transformer(self.paw.density.finegd, wfs.gd,
                                    self.paw.input_parameters.stencils[1]
                                    ).apply

        if xc == 'RPA':
            xc = None  # enable RPA as keyword
        if xc is not None:
            self.xc = XC(xc)
            self.xc.initialize(self.paw.density, self.paw.hamiltonian,
                               wfs, self.paw.occupations)

            # check derivativeLevel
            if derivativeLevel is None:
                derivativeLevel = \
                    self.xc.get_functional().get_max_derivative_level()
            self.derivativeLevel = derivativeLevel
            # change the setup xc functional if needed
            # the ground state calculation may have used another xc
            if kss.npspins > kss.nvspins:
                spin_increased = True
            else:
                spin_increased = False
        else:
            self.xc = None

        self.numscale = numscale

        self.singletsinglet = False
        if kss.nvspins < 2 and kss.npspins < 2:
            # this will be a singlet to singlet calculation only
            self.singletsinglet = True

        nij = len(kss)
        self.Om = np.zeros((nij, nij))
        self.get_full()

    def get_full(self):

        self.paw.timer.start('Omega RPA')
        self.get_rpa()
        self.paw.timer.stop()

        if self.xc is not None:
            self.paw.timer.start('Omega XC')
            self.get_xc()
            self.paw.timer.stop()

        self.eh_comm.sum(self.Om)
        self.full = self.Om

    def get_xc(self):
        """Add xc part of the coupling matrix"""

        # shorthands
        paw = self.paw
        wfs = paw.wfs
        gd = paw.density.finegd
        comm = gd.comm
        eh_comm = self.eh_comm

        fg = self.finegrid is 2
        kss = self.fullkss
        nij = len(kss)

        Om_xc = self.Om
        # initialize densities
        # nt_sg is the smooth density on the fine grid with spin index

        if kss.nvspins == 2:
            # spin polarised ground state calc.
            nt_sg = paw.density.nt_sg
        else:
            # spin unpolarised ground state calc.
            if kss.npspins == 2:
                # construct spin polarised densities
                nt_sg = np.array([.5 * paw.density.nt_sg[0],
                                  .5 * paw.density.nt_sg[0]])
            else:
                nt_sg = paw.density.nt_sg
        # check if D_sp have been changed before
        D_asp = self.paw.density.D_asp
        for a, D_sp in D_asp.items():
            if len(D_sp) != kss.npspins:
                if len(D_sp) == 1:
                    D_asp[a] = np.array([0.5 * D_sp[0], 0.5 * D_sp[0]])
                else:
                    D_asp[a] = np.array([D_sp[0] + D_sp[1]])

        # restrict the density if needed
        if fg:
            nt_s = nt_sg
        else:
            nt_s = self.gd.zeros(nt_sg.shape[0])
            for s in range(nt_sg.shape[0]):
                self.restrict(nt_sg[s], nt_s[s])
            gd = paw.density.gd

        # initialize vxc or fxc

        if self.derivativeLevel == 0:
            raise NotImplementedError
            if kss.npspins == 2:
                v_g = nt_sg[0].copy()
            else:
                v_g = nt_sg.copy()
        elif self.derivativeLevel == 1:
            pass
        elif self.derivativeLevel == 2:
            fxc_sg = np.zeros(nt_sg.shape)
            self.xc.calculate_fxc(gd, nt_sg, fxc_sg)
        else:
            raise ValueError('derivativeLevel can only be 0,1,2')

# self.paw.my_nuclei = []

        ns = self.numscale
        xc = self.xc
        print('XC', nij, 'transitions', file=self.txt)
        for ij in range(eh_comm.rank, nij, eh_comm.size):
            print('XC kss[' + '%d' % ij + ']', file=self.txt)

            timer = Timer()
            timer.start('init')
            timer2 = Timer()

            if self.derivativeLevel >= 1:
                # vxc is available
                # We use the numerical two point formula for calculating
                # the integral over fxc*n_ij. The results are
                # vvt_s        smooth integral
                # nucleus.I_sp atom based correction matrices (pack2)
                #              stored on each nucleus
                timer2.start('init v grids')
                vp_s = np.zeros(nt_s.shape, nt_s.dtype.char)
                vm_s = np.zeros(nt_s.shape, nt_s.dtype.char)
                if kss.npspins == 2:  # spin polarised
                    nv_s = nt_s.copy()
                    nv_s[kss[ij].pspin] += ns * kss[ij].get(fg)
                    xc.calculate(gd, nv_s, vp_s)
                    nv_s = nt_s.copy()
                    nv_s[kss[ij].pspin] -= ns * kss[ij].get(fg)
                    xc.calculate(gd, nv_s, vm_s)
                else:  # spin unpolarised
                    nv = nt_s + ns * kss[ij].get(fg)
                    xc.calculate(gd, nv, vp_s)
                    nv = nt_s - ns * kss[ij].get(fg)
                    xc.calculate(gd, nv, vm_s)
                vvt_s = (0.5 / ns) * (vp_s - vm_s)
                timer2.stop()

                # initialize the correction matrices
                timer2.start('init v corrections')
                I_asp = {}
                for a, P_ni in wfs.kpt_u[kss[ij].spin].P_ani.items():
                    # create the modified density matrix
                    Pi_i = P_ni[kss[ij].i]
                    Pj_i = P_ni[kss[ij].j]
                    P_ii = np.outer(Pi_i, Pj_i)
                    # we need the symmetric form, hence we can pack
                    P_p = pack(P_ii)
                    D_sp = self.paw.density.D_asp[a].copy()
                    D_sp[kss[ij].pspin] -= ns * P_p
                    setup = wfs.setups[a]
                    I_sp = np.zeros_like(D_sp)
                    self.xc.calculate_paw_correction(setup, D_sp, I_sp)
                    I_sp *= -1.0
                    D_sp = self.paw.density.D_asp[a].copy()
                    D_sp[kss[ij].pspin] += ns * P_p
                    self.xc.calculate_paw_correction(setup, D_sp, I_sp)
                    I_sp /= 2.0 * ns
                    I_asp[a] = I_sp
                timer2.stop()

            timer.stop()
            t0 = timer.get_time('init')
            timer.start(ij)

            for kq in range(ij, nij):
                weight = self.weight_Kijkq(ij, kq)

                if self.derivativeLevel == 0:
                    # only Exc is available

                    if kss.npspins == 2:  # spin polarised
                        nv_g = nt_sg.copy()
                        nv_g[kss[ij].pspin] += kss[ij].get(fg)
                        nv_g[kss[kq].pspin] += kss[kq].get(fg)
                        Excpp = xc.get_energy_and_potential(
                            nv_g[0], v_g, nv_g[1], v_g)
                        nv_g = nt_sg.copy()
                        nv_g[kss[ij].pspin] += kss[ij].get(fg)
                        nv_g[kss[kq].pspin] -= kss[kq].get(fg)
                        Excpm = xc.get_energy_and_potential(
                            nv_g[0], v_g, nv_g[1], v_g)
                        nv_g = nt_sg.copy()
                        nv_g[kss[ij].pspin] -=\
                            kss[ij].get(fg)
                        nv_g[kss[kq].pspin] +=\
                            kss[kq].get(fg)
                        Excmp = xc.get_energy_and_potential(
                            nv_g[0], v_g, nv_g[1], v_g)
                        nv_g = nt_sg.copy()
                        nv_g[kss[ij].pspin] -= \
                            kss[ij].get(fg)
                        nv_g[kss[kq].pspin] -=\
                            kss[kq].get(fg)
                        Excpp = xc.get_energy_and_potential(
                            nv_g[0], v_g, nv_g[1], v_g)
                    else:  # spin unpolarised
                        nv_g = nt_sg + ns * kss[ij].get(fg)\
                            + ns * kss[kq].get(fg)
                        Excpp = xc.get_energy_and_potential(nv_g, v_g)
                        nv_g = nt_sg + ns * kss[ij].get(fg)\
                            - ns * kss[kq].get(fg)
                        Excpm = xc.get_energy_and_potential(nv_g, v_g)
                        nv_g = nt_sg - ns * kss[ij].get(fg)\
                            + ns * kss[kq].get(fg)
                        Excmp = xc.get_energy_and_potential(nv_g, v_g)
                        nv_g = nt_sg - ns * kss[ij].get(fg)\
                            - ns * kss[kq].get(fg)
                        Excmm = xc.get_energy_and_potential(nv_g, v_g)

                    Om_xc[ij, kq] += weight *\
                        0.25 * \
                        (Excpp - Excmp - Excpm + Excmm) / (ns * ns)

                elif self.derivativeLevel == 1:
                    # vxc is available

                    timer2.start('integrate')
                    Om_xc[ij, kq] += weight *\
                        self.gd.integrate(kss[kq].get(fg) *
                                          vvt_s[kss[kq].pspin])
                    timer2.stop()

                    timer2.start('integrate corrections')
                    Exc = 0.
                    for a, P_ni in wfs.kpt_u[kss[kq].spin].P_ani.items():
                        # create the modified density matrix
                        Pk_i = P_ni[kss[kq].i]
                        Pq_i = P_ni[kss[kq].j]
                        P_ii = np.outer(Pk_i, Pq_i)
                        # we need the symmetric form, hence we can pack
                        # use pack as I_sp used pack2
                        P_p = pack(P_ii)
                        Exc += np.dot(I_asp[a][kss[kq].pspin], P_p)
                    Om_xc[ij, kq] += weight * self.gd.comm.sum(Exc)
                    timer2.stop()

                elif self.derivativeLevel == 2:
                    # fxc is available
                    if kss.npspins == 2:  # spin polarised
                        Om_xc[ij, kq] += weight *\
                            gd.integrate(kss[ij].get(fg) *
                                         kss[kq].get(fg) *
                                         fxc_sg[kss[ij].pspin, kss[kq].pspin])
                    else:  # spin unpolarised
                        Om_xc[ij, kq] += weight *\
                            gd.integrate(kss[ij].get(fg) *
                                         kss[kq].get(fg) *
                                         fxc_sg)

                    # XXX still numeric derivatives for local terms
                    timer2.start('integrate corrections')
                    Exc = 0.
                    for a, P_ni in wfs.kpt_u[kss[kq].spin].P_ani.items():
                        # create the modified density matrix
                        Pk_i = P_ni[kss[kq].i]
                        Pq_i = P_ni[kss[kq].j]
                        P_ii = np.outer(Pk_i, Pq_i)
                        # we need the symmetric form, hence we can pack
                        # use pack as I_sp used pack2
                        P_p = pack(P_ii)
                        Exc += np.dot(I_asp[a][kss[kq].pspin], P_p)
                    Om_xc[ij, kq] += weight * self.gd.comm.sum(Exc)
                    timer2.stop()

                if ij != kq:
                    Om_xc[kq, ij] = Om_xc[ij, kq]

            timer.stop()
# timer2.write()
            if ij < (nij - 1):
                print('XC estimated time left',
                      self.time_left(timer, t0, ij, nij), file=self.txt)

    def Coulomb_integral_kss(self, kss_ij, kss_kq, phit, rhot,
                             timer=None):
        # smooth part
        if timer:
            timer.start('integrate')
        I = self.gd.integrate(rhot * phit)
        if timer:
            timer.stop()
            timer.start('integrate corrections 2')

        wfs = self.paw.wfs
        Pij_ani = wfs.kpt_u[kss_ij.spin].P_ani
        Pkq_ani = wfs.kpt_u[kss_kq.spin].P_ani

        # Add atomic corrections
        Ia = 0.0
        for a, Pij_ni in Pij_ani.items():
            Pi_i = Pij_ni[kss_ij.i]
            Pj_i = Pij_ni[kss_ij.j]
            Dij_ii = np.outer(Pi_i, Pj_i)
            Dij_p = pack(Dij_ii)
            Pk_i = Pkq_ani[a][kss_kq.i]
            Pq_i = Pkq_ani[a][kss_kq.j]
            Dkq_ii = np.outer(Pk_i, Pq_i)
            Dkq_p = pack(Dkq_ii)
            C_pp = wfs.setups[a].M_pp
            #   ----
            # 2 >      P   P  C    P  P
            #   ----    ip  jr prst ks qt
            #   prst
            Ia += 2.0 * np.dot(Dkq_p, np.dot(C_pp, Dij_p))
        I += self.gd.comm.sum(Ia)
        if timer:
            timer.stop()

        return I

    def get_rpa(self):
        """calculate RPA part of the omega matrix"""

        # shorthands
        kss = self.fullkss
        finegrid = self.finegrid
        wfs = self.paw.wfs
        eh_comm = self.eh_comm

        # calculate omega matrix
        nij = len(kss)
        print('RPA', nij, 'transitions', file=self.txt)

        Om = self.Om

        for ij in range(eh_comm.rank, nij, eh_comm.size):
            print('RPA kss[' + '%d' % ij + ']=', kss[ij], file=self.txt)

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
            phit_p = np.zeros(rhot_p.shape, rhot_p.dtype.char)
            self.poisson.solve(phit_p, rhot_p, charge=None)
            timer2.stop()

            timer.stop()
            t0 = timer.get_time('init')
            timer.start(ij)

            if finegrid == 1:
                rhot = kss[ij].with_compensation_charges()
                phit = self.gd.zeros()
# print "shapes 0=",phit.shape,rhot.shape
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

                pre = 2 * sqrt(kss[ij].get_energy() * kss[kq].get_energy() *
                               kss[ij].get_weight() * kss[kq].get_weight())
                I = self.Coulomb_integral_kss(kss[ij], kss[kq],
                                              rhot, phit, timer2)
                Om[ij, kq] = pre * I

                if ij == kq:
                    Om[ij, kq] += kss[ij].get_energy() ** 2
                else:
                    Om[kq, ij] = Om[ij, kq]

            timer.stop()
# timer2.write()
            if ij < (nij - 1):
                t = timer.get_time(ij)  # time for nij-ij calculations
                t = .5 * t * \
                    (nij - ij)  # estimated time for n*(n+1)/2, n=nij-(ij+1)
                print('RPA estimated time left',
                      self.timestring(t0 * (nij - ij - 1) + t), file=self.txt)

    def singlets_triplets(self):
        """Split yourself into singlet and triplet transitions"""

        assert(self.fullkss.npspins == 2)
        assert(self.fullkss.nvspins == 1)

        # strip kss from down spins
        skss = KSSingles()
        tkss = KSSingles()
        map = []
        for ij, ks in enumerate(self.fullkss):
            if ks.pspin == ks.spin:
                skss.append((ks + ks) / sqrt(2))
                tkss.append((ks - ks) / sqrt(2))
                map.append(ij)

        nkss = len(skss)

        # define the singlet and the triplet omega-matrixes
        sOm = OmegaMatrix(kss=skss)
        sOm.full = np.empty((nkss, nkss))
        tOm = OmegaMatrix(kss=tkss)
        tOm.full = np.empty((nkss, nkss))
        for ij in range(nkss):
            for kl in range(nkss):
                sOm.full[ij, kl] = (self.full[map[ij], map[kl]] +
                                    self.full[map[ij], nkss + map[kl]])
                tOm.full[ij, kl] = (self.full[map[ij], map[kl]] -
                                    self.full[map[ij], nkss + map[kl]])
        return sOm, tOm

    def timestring(self, t):
        ti = int(t + 0.5)
        td = ti // 86400
        st = ''
        if td > 0:
            st += '%d' % td + 'd'
            ti -= td * 86400
        th = ti // 3600
        if th > 0:
            st += '%d' % th + 'h'
            ti -= th * 3600
        tm = ti // 60
        if tm > 0:
            st += '%d' % tm + 'm'
            ti -= tm * 60
        st += '%d' % ti + 's'
        return st

    def time_left(self, timer, t0, ij, nij):
        t = timer.get_time(ij)  # time for nij-ij calculations
        t = .5 * t * (nij - ij)  # estimated time for n*(n+1)/2, n=nij-(ij+1)
        return self.timestring(t0 * (nij - ij - 1) + t)

    def get_map(self, istart=None, jend=None, energy_range=None):
        """Return the reduction map for the given requirements"""

        self.istart = istart
        self.jend = jend
        if istart is None and jend is None and energy_range is None:
            return None, self.fullkss

        # reduce the matrix
        print('# diagonalize: %d transitions original'
              % len(self.fullkss), file=self.txt)

        if energy_range is None:
            if istart is None:
                istart = self.kss.istart
            if self.fullkss.istart > istart:
                raise RuntimeError('istart=%d has to be >= %d' %
                                   (istart, self.kss.istart))
            if jend is None:
                jend = self.kss.jend
            if self.fullkss.jend < jend:
                raise RuntimeError('jend=%d has to be <= %d' %
                                   (jend, self.kss.jend))

        else:
            try:
                emin, emax = energy_range
            except:
                emax = energy_range
                emin = 0.
            emin /= Hartree
            emax /= Hartree

        map = []
        kss = KSSingles()
        for ij, k in zip(range(len(self.fullkss)), self.fullkss):
            if energy_range is None:
                if k.i >= istart and k.j <= jend:
                    kss.append(k)
                    map.append(ij)
            else:
                if k.energy >= emin and k.energy < emax:
                    kss.append(k)
                    map.append(ij)
        kss.update()
        print('# diagonalize: %d transitions now' % len(kss), file=self.txt)

        return map, kss

    def diagonalize(self, istart=None, jend=None, energy_range=None,
                    TDA=False):
        """Evaluate Eigenvectors and Eigenvalues:"""

        if TDA:
            raise NotImplementedError

        map, kss = self.get_map(istart, jend, energy_range)
        nij = len(kss)
        if map is None:
            evec = self.full.copy()
        else:
            evec = np.zeros((nij, nij))
            for ij in range(nij):
                for kq in range(nij):
                    evec[ij, kq] = self.full[map[ij], map[kq]]
        assert(len(evec) > 0)

        self.eigenvectors = evec
        self.eigenvalues = np.zeros((len(kss)))
        self.kss = kss
        diagonalize(self.eigenvectors, self.eigenvalues)

    def Kss(self, kss=None):
        """Set and get own Kohn-Sham singles"""
        if kss is not None:
            self.fullkss = kss
        if(hasattr(self, 'fullkss')):
            return self.fullkss
        else:
            return None

    def read(self, filename=None, fh=None):
        """Read myself from a file"""
        if fh is None:
            f = open(filename, 'r')
        else:
            f = fh

        f.readline()
        nij = int(f.readline())
        full = np.zeros((nij, nij))
        for ij in range(nij):
            l = f.readline().split()
            for kq in range(ij, nij):
                full[ij, kq] = float(l[kq - ij])
                full[kq, ij] = full[ij, kq]
        self.full = full

        if fh is None:
            f.close()

    def write(self, filename=None, fh=None):
        """Write current state to a file."""
        if mpi.rank == mpi.MASTER:
            if fh is None:
                f = open(filename, 'w')
            else:
                f = fh

            f.write('# OmegaMatrix\n')
            nij = len(self.fullkss)
            f.write('%d\n' % nij)
            for ij in range(nij):
                for kq in range(ij, nij):
                    f.write(' %g' % self.full[ij, kq])
                f.write('\n')

            if fh is None:
                f.close()

    def weight_Kijkq(self, ij, kq):
        """weight for the coupling matrix terms"""
        kss = self.fullkss
        return 2. * sqrt(kss[ij].get_energy() * kss[kq].get_energy() *
                         kss[ij].get_weight() * kss[kq].get_weight())

    def __str__(self):
        str = '<OmegaMatrix> '
        if hasattr(self, 'eigenvalues'):
            str += 'dimension ' + ('%d' % len(self.eigenvalues))
            str += '\neigenvalues: '
            for ev in self.eigenvalues:
                str += ' ' + ('%f' % (sqrt(ev) * Hartree))
        return str
