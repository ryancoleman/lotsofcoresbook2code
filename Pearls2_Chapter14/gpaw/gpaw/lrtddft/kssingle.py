"""Kohn-Sham single particle excitations realated objects.

"""
from __future__ import print_function
import sys
from math import pi, sqrt

import numpy as np
from ase.units import Bohr, Hartree, alpha
from ase.parallel import paropen

import gpaw.mpi as mpi
from gpaw.utilities import packed_index
from gpaw.lrtddft.excitation import Excitation, ExcitationList
from gpaw.pair_density import PairDensity
from gpaw.fd_operators import Gradient
from gpaw.utilities.tools import coordinates


class KSSingles(ExcitationList):

    """Kohn-Sham single particle excitations

    Input parameters:

    calculator:
      the calculator object after a ground state calculation

    nspins:
      number of spins considered in the calculation
      Note: Valid only for unpolarised ground state calculation

    eps:
      Minimal occupation difference for a transition (default 0.001)

    istart:
      First occupied state to consider
    jend:
      Last unoccupied state to consider
    energy_range:
      The energy range [emin, emax] or emax for KS transitions to use as basis
    """

    def __init__(self,
                 calculator=None,
                 nspins=None,
                 eps=0.001,
                 istart=0,
                 jend=sys.maxsize,
                 energy_range=None,
                 filehandle=None,
                 txt=None):

        self.eps = None

        if isinstance(calculator, str):
            filehandle = open(calculator)
        if filehandle is not None:
            self.world = mpi.world
            self.read(fh=filehandle, istart=istart, jend=jend)
            return None

        # LCAO calculation requires special actions
        if calculator is not None:
            self.lcao = calculator.input_parameters.mode == 'lcao'

        ExcitationList.__init__(self, calculator, txt=txt)

        if calculator is None:
            return  # leave the list empty

        # deny hybrids as their empty states are wrong
#        gsxc = calculator.hamiltonian.xc
#        hybrid = hasattr(gsxc, 'hybrid') and gsxc.hybrid > 0.0
#        assert(not hybrid)

        # ensure correctly initialized wave functions
        calculator.converge_wave_functions()
        self.world = calculator.wfs.world

        # parallelization over bands not yet supported
        assert(calculator.wfs.band_comm.size == 1)

        self.select(nspins, eps, istart, jend, energy_range)

        trkm = self.get_trk()
        print('KSS TRK sum %g (%g,%g,%g)' %
              (np.sum(trkm) / 3., trkm[0], trkm[1], trkm[2]), file=self.txt)
        pol = self.get_polarizabilities(lmax=3)
        print('KSS polarisabilities(l=0-3) %g, %g, %g, %g' %
              tuple(pol.tolist()), file=self.txt)

    def select(self, nspins=None, eps=0.001,
               istart=0, jend=sys.maxsize, energy_range=None):
        """Select KSSingles according to the given criterium."""

        paw = self.calculator
        wfs = paw.wfs
        self.dtype = wfs.dtype
        self.kpt_u = wfs.kpt_u

        if not self.lcao and self.kpt_u[0].psit_nG is None:
            raise RuntimeError('No wave functions in calculator!')

        # criteria
        emin = -sys.float_info.max
        emax = sys.float_info.max
        if energy_range is not None:
            try:
                emin, emax = energy_range
                emin /= Hartree
                emax /= Hartree
            except:
                emax = energy_range / Hartree
        self.istart = istart
        self.jend = jend
        self.eps = eps

        # here, we need to take care of the spins also for
        # closed shell systems (Sz=0)
        # vspin is the virtual spin of the wave functions,
        #       i.e. the spin used in the ground state calculation
        # pspin is the physical spin of the wave functions
        #       i.e. the spin of the excited states
        self.nvspins = wfs.nspins
        self.npspins = wfs.nspins
        fijscale = 1
        ispins = [0]
        nks = wfs.kd.nks
        if self.nvspins < 2:
            if nspins > self.nvspins:
                self.npspins = nspins
                fijscale = 0.5
                ispins = [0, 1]
                nks = 2 * wfs.kd.nks

        kpt_comm = self.calculator.wfs.kd.comm
        nbands = len(self.kpt_u[0].f_n)

        # select
        take = np.zeros((nks, nbands, nbands), dtype=int)
        u = 0
        for ispin in ispins:
            for ks in range(wfs.kd.nks):
                myks = ks - wfs.kd.ks0
                if myks >= 0 and myks < wfs.kd.mynks:
                    kpt = self.kpt_u[myks]
                    for i in range(nbands):
                        for j in range(i + 1, nbands):
                            fij = kpt.f_n[i] - kpt.f_n[j]
                            epsij = kpt.eps_n[j] - kpt.eps_n[i]
                            if (fij > eps and
                                epsij >= emin and epsij < emax and
                                    i >= self.istart and j <= self.jend):
                                take[u, i, j] = 1
                u += 1
        kpt_comm.sum(take)

        # calculate in parallel
        u = 0
        for ispin in ispins:
            for ks in range(wfs.kd.nks):
                myks = ks - wfs.kd.ks0
                for i in range(nbands):
                    for j in range(i + 1, nbands):
                        if take[u, i, j]:
                            if myks >= 0 and myks < wfs.kd.mynks:
                                kpt = self.kpt_u[myks]
                                pspin = max(kpt.s, ispin)
                                self.append(
                                    KSSingle(i, j, pspin, kpt, paw,
                                             fijscale=fijscale,
                                             dtype=self.dtype))
                            else:
                                self.append(KSSingle(i, j, pspin=0,
                                                     kpt=None, paw=paw,
                                                     dtype=self.dtype))
                u += 1

        # distribute
        for kss in self:
            kss.distribute()

    def read(self, filename=None, fh=None, istart=0, jend=sys.maxsize):
        """Read myself from a file"""
        if fh is None:
            if filename.endswith('.gz'):
                import gzip
                f = gzip.open(filename)
            else:
                f = open(filename, 'r')
        else:
            f = fh

        try:
            assert(f.readline().strip() == '# KSSingles')
        except:
            raise RuntimeError(f.name + ' is not a ' +
                               self.__class__.__name__ + ' data file')
        words = f.readline().split()
        n = int(words[0])
        if len(words) == 1:
            # old output style for real wave functions (finite systems)
            self.dtype = float
        else:
            if words[1].startswith('complex'):
                self.dtype = complex
            else:
                self.dtype = float
            self.eps = float(f.readline())
        self.npspins = 1
        for i in range(n):
            kss = KSSingle(string=f.readline(), dtype=self.dtype)
            if (kss.i >= istart) and (kss.j <= jend):
                self.append(kss)
                self.npspins = max(self.npspins, kss.pspin + 1)
        self.update()

        if fh is None:
            f.close()

    def update(self):
        istart = self[0].i
        jend = 0
        npspins = 1
        nvspins = 1
        for kss in self:
            istart = min(kss.i, istart)
            jend = max(kss.j, jend)
            if kss.pspin == 1:
                npspins = 2
            if kss.spin == 1:
                nvspins = 2
        self.istart = istart
        self.jend = jend
        self.npspins = npspins
        self.nvspins = nvspins

        if hasattr(self, 'energies'):
            del(self.energies)

    def set_arrays(self):
        if hasattr(self, 'energies'):
            return
        energies = []
        fij = []
        me = []
        mur = []
        muv = []
        magn = []
        for k in self:
            energies.append(k.energy)
            fij.append(k.fij)
            me.append(k.me)
            mur.append(k.mur)
            if k.muv is not None:
                muv.append(k.muv)
            if k.magn is not None:
                magn.append(k.magn)
        self.energies = np.array(energies)
        self.fij = np.array(fij)
        self.me = np.array(me)
        self.mur = np.array(mur)
        if len(muv):
            self.muv = np.array(muv)
        else:
            self.muv = None
        if len(magn):
            self.magn = np.array(magn)
        else:
            self.magn = None

    def write(self, filename=None, fh=None):
        """Write current state to a file.

        'filename' is the filename. If the filename ends in .gz,
        the file is automatically saved in compressed gzip format.

        'fh' is a filehandle. This can be used to write into already
        opened files.
        """
        if self.world.rank != 0:
            return

        if fh is None:
            if filename.endswith('.gz') and mpi.rank == mpi.MASTER:
                import gzip
                f = gzip.open(filename, 'wb')
            else:
                f = paropen(filename, 'w')
        else:
            f = fh

        f.write('# KSSingles\n')
        f.write('{0} {1}\n'.format(len(self), np.dtype(self.dtype)))
        f.write('{0}\n'.format(self.eps))
        for kss in self:
            f.write(kss.outstring())

        if fh is None:
            f.close()


class KSSingle(Excitation, PairDensity):

    """Single Kohn-Sham transition containing all it's indicees

    ::

      pspin=physical spin
      spin=virtual  spin, i.e. spin in the ground state calc.
      kpt=the Kpoint object
      fijscale=weight for the occupation difference::
      me  = sqrt(fij*epsij) * <i|r|j>
      mur = - <i|r|a>
      muv = - <i|nabla|a>/omega_ia with omega_ia>0
      magn = <i|[r x nabla]|a> / (2 m_e c)
    """

    def __init__(self, iidx=None, jidx=None, pspin=None, kpt=None,
                 paw=None, string=None, fijscale=1, dtype=float):

        if string is not None:
            self.fromstring(string, dtype)
            return None

        # normal entry

        PairDensity.__init__(self, paw)
        PairDensity.initialize(self, kpt, iidx, jidx)

        self.pspin = pspin

        self.energy = 0.0
        self.fij = 0.0

        self.me = np.zeros((3), dtype=dtype)
        self.mur = np.zeros((3), dtype=dtype)
        self.muv = np.zeros((3), dtype=dtype)
        self.magn = np.zeros((3), dtype=dtype)

        self.kpt_comm = paw.wfs.kd.comm

        # leave empty if not my kpt
        if kpt is None:
            return

        wfs = paw.wfs
        gd = wfs.gd

        self.energy = kpt.eps_n[jidx] - kpt.eps_n[iidx]
        self.fij = (kpt.f_n[iidx] - kpt.f_n[jidx]) * fijscale

        # calculate matrix elements -----------

        # length form ..........................

        # course grid contribution
        # <i|r|j> is the negative of the dipole moment (because of negative
        # e- charge)
        me = - gd.calculate_dipole_moment(self.get())

        # augmentation contributions
        ma = np.zeros(me.shape, dtype=dtype)
        pos_av = paw.atoms.get_positions() / Bohr
        for a, P_ni in kpt.P_ani.items():
            Ra = pos_av[a]
            Pi_i = P_ni[self.i].conj()
            Pj_i = P_ni[self.j]
            Delta_pL = wfs.setups[a].Delta_pL
            ni = len(Pi_i)
            ma0 = 0
            ma1 = np.zeros(me.shape, dtype=me.dtype)
            for i in range(ni):
                for j in range(ni):
                    pij = Pi_i[i] * Pj_i[j]
                    ij = packed_index(i, j, ni)
                    # L=0 term
                    ma0 += Delta_pL[ij, 0] * pij
                    # L=1 terms
                    if wfs.setups[a].lmax >= 1:
                        # see spherical_harmonics.py for
                        # L=1:y L=2:z; L=3:x
                        ma1 += np.array([Delta_pL[ij, 3], Delta_pL[ij, 1],
                                         Delta_pL[ij, 2]]) * pij
            ma += sqrt(4 * pi / 3) * ma1 + Ra * sqrt(4 * pi) * ma0
        gd.comm.sum(ma)

        self.me = sqrt(self.energy * self.fij) * (me + ma)
        self.mur = - (me + ma)

        # velocity form .............................

        if self.lcao:
            # XXX Velocity form not supported in LCAO
            return

        me = np.zeros(self.mur.shape, dtype=dtype)

        # get derivatives
        dtype = self.wfj.dtype
        dwfj_cg = gd.empty((3), dtype=dtype)
        if not hasattr(gd, 'ddr'):
            gd.ddr = [Gradient(gd, c, dtype=dtype).apply for c in range(3)]
        for c in range(3):
            gd.ddr[c](self.wfj, dwfj_cg[c], kpt.phase_cd)
            me[c] = gd.integrate(self.wfi.conj() * dwfj_cg[c])

        if 0:
            me2 = np.zeros(self.mur.shape)
            for c in range(3):
                gd.ddr[c](self.wfi, dwfj_cg[c], kpt.phase_cd)
                me2[c] = gd.integrate(self.wfj * dwfj_cg[c])
            print(me, -me2, me2 + me)

        # augmentation contributions
        ma = np.zeros(me.shape, dtype=me.dtype)
        for a, P_ni in kpt.P_ani.items():
            Pi_i = P_ni[self.i].conj()
            Pj_i = P_ni[self.j]
            nabla_iiv = paw.wfs.setups[a].nabla_iiv
            for c in range(3):
                for i1, Pi in enumerate(Pi_i):
                    for i2, Pj in enumerate(Pj_i):
                        ma[c] += Pi * Pj * nabla_iiv[i1, i2, c]
        gd.comm.sum(ma)

        self.muv = - (me + ma) / self.energy

        # magnetic transition dipole ................

        r_cg, r2_g = coordinates(gd)
        magn = np.zeros(me.shape, dtype=dtype)

        wfi_g = self.wfi.conj()
        for ci in range(3):
            cj = (ci + 1) % 3
            ck = (ci + 2) % 3
            magn[ci] = gd.integrate(wfi_g * r_cg[cj] * dwfj_cg[ck] -
                                    wfi_g * r_cg[ck] * dwfj_cg[cj])
        # augmentation contributions
        ma = np.zeros(magn.shape, dtype=magn.dtype)
        for a, P_ni in kpt.P_ani.items():
            Pi_i = P_ni[self.i].conj()
            Pj_i = P_ni[self.j]
            rnabla_iiv = paw.wfs.setups[a].rnabla_iiv
            for c in range(3):
                for i1, Pi in enumerate(Pi_i):
                    for i2, Pj in enumerate(Pj_i):
                        ma[c] += Pi * Pj * rnabla_iiv[i1, i2, c]
        gd.comm.sum(ma)

        self.magn = -alpha / 2. * (magn + ma)

    def distribute(self):
        """Distribute results to all cores."""
        self.spin = self.kpt_comm.sum(self.spin)
        self.pspin = self.kpt_comm.sum(self.pspin)
        self.k = self.kpt_comm.sum(self.k)
        self.weight = self.kpt_comm.sum(self.weight)
        self.energy = self.kpt_comm.sum(self.energy)
        self.fij = self.kpt_comm.sum(self.fij)

        self.kpt_comm.sum(self.me)
        self.kpt_comm.sum(self.mur)
        self.kpt_comm.sum(self.muv)
        self.kpt_comm.sum(self.magn)

    def __add__(self, other):
        """Add two KSSingles"""
        result = self.copy()
        result.me = self.me + other.me
        result.mur = self.mur + other.mur
        result.muv = self.muv + other.muv
        result.magn = self.magn - other.magn
        return result

    def __sub__(self, other):
        """Subtract two KSSingles"""
        result = self.copy()
        result.me = self.me - other.me
        result.mur = self.mur - other.mur
        result.muv = self.muv - other.muv
        result.magn = self.magn - other.magn
        return result

    def __rmul__(self, x):
        return self.__mul__(x)

    def __mul__(self, x):
        """Multiply a KSSingle with a number"""
        if isinstance(x, (float, int)):
            result = self.copy()
            result.me = self.me * x
            result.mur = self.mur * x
            result.muv = self.muv * x
            return result
        else:
            return RuntimeError('not a number')

    def __div__(self, x):
        return self.__mul__(1. / x)

    def copy(self):
        if self.mur.dtype == complex:
            return KSSingle(string=self.outstring(), dtype=complex)
        else:
            return KSSingle(string=self.outstring(), dtype=float)

    def fromstring(self, string, dtype=float):
        l = string.split()
        self.i = int(l.pop(0))
        self.j = int(l.pop(0))
        self.pspin = int(l.pop(0))
        self.spin = int(l.pop(0))
        if dtype == float:
            self.k = 0
            self.weight = 1
        else:
            self.k = int(l.pop(0))
            self.weight = float(l.pop(0))
        self.energy = float(l.pop(0))
        self.fij = float(l.pop(0))
        self.mur = np.array([dtype(l.pop(0)) for i in range(3)])
        self.me = - self.mur * sqrt(self.energy * self.fij)
        self.muv = self.magn = None
        if len(l):
            self.muv = np.array([dtype(l.pop(0)) for i in range(3)])
        if len(l):
            self.magn = np.array([dtype(l.pop(0)) for i in range(3)])
        return None

    def outstring(self):
        if self.mur.dtype == float:
            string = '{0:d} {1:d}  {2:d} {3:d}  {4:g} {5:g}'.format(
                self.i, self.j, self.pspin, self.spin, self.energy, self.fij)
        else:
            string = (
                '{0:d} {1:d}  {2:d} {3:d} {4:d} {5:g}  {6:g} {7:g}'.format(
                    self.i, self.j, self.pspin, self.spin, self.k,
                    self.weight, self.energy, self.fij))
        string += '  '

        def format_me(me):
            string = ''
            if me.dtype == float:
                for m in me:
                    string += ' {0:.5e}'.format(m)
            else:
                for m in me:
                    string += ' {0.real:.5e}{0.imag:+.5e}j'.format(m)
            return string

        string += '  ' + format_me(self.mur)
        if self.muv is not None:
            string += '  ' + format_me(self.muv)
        if self.magn is not None:
            string += '  ' + format_me(self.magn)
        string += '\n'

        return string

    def __str__(self):
        string = '# <KSSingle> %d->%d %d(%d) eji=%g[eV]' % \
            (self.i, self.j, self.pspin, self.spin,
             self.energy * Hartree)
        if self.me.dtype == float:
            string += ' (%g,%g,%g)' % (self.me[0], self.me[1], self.me[2])
        else:
            string += ' kpt={0:d} w={1:g}'.format(self.k, self.weight)
            string += ' ('
            # use velocity form
            s = - np.sqrt(self.energy * self.fij)
            for c, m in enumerate(s * self.me):
                string += '{0.real:.5e}{0.imag:+.5e}j'.format(m)
                if c < 2:
                    string += ','
            string += ')'
        return string

    #
    # User interface: ##
    #

    def get_weight(self):
        return self.fij
