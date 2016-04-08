"""This module defines a linear response TDDFT-class.

"""
from __future__ import print_function
from math import sqrt
import sys

import numpy as np
from ase.units import Hartree
from ase.utils.timing import Timer

import _gpaw
import gpaw.mpi as mpi
MASTER = mpi.MASTER
from gpaw import debug
from gpaw.lrtddft.excitation import Excitation, ExcitationList
from gpaw.lrtddft.kssingle import KSSingles
from gpaw.lrtddft.omega_matrix import OmegaMatrix
from gpaw.lrtddft.apmb import ApmB
# from gpaw.lrtddft.transition_density import TransitionDensity
from gpaw.utilities import packed_index
from gpaw.utilities.lapack import diagonalize
from gpaw.xc import XC
from gpaw.xc.hybridk import HybridXC
from gpaw.utilities.timing import nulltimer
from gpaw.lrtddft.spectrum import spectrum
from gpaw.wavefunctions.fd import FDWaveFunctions

__all__ = ['LrTDDFT', 'photoabsorption_spectrum', 'spectrum']


class LrTDDFT(ExcitationList):

    """Linear Response TDDFT excitation class

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

    xc:
    Exchange-Correlation approximation in the Kernel
    derivative_level:
    0: use Exc, 1: use vxc, 2: use fxc  if available

    filename:
    read from a file
    """

    def __init__(self, calculator=None, **kwargs):

        self.timer = Timer()

        self.set(**kwargs)

        if isinstance(calculator, str):
            ExcitationList.__init__(self, None, self.txt)
            self.filename = calculator
        else:
            ExcitationList.__init__(self, calculator, self.txt)

        if self.filename is not None:
            return self.read(self.filename)

        if self.eh_comm is None:
            self.eh_comm = mpi.serial_comm
        elif isinstance(self.eh_comm, (mpi.world.__class__,
                                       mpi.serial_comm.__class__)):
            # Correct type already.
            pass
        else:
            # world should be a list of ranks:
            self.eh_comm = mpi.world.new_communicator(np.asarray(eh_comm))

        if calculator is not None and calculator.initialized:
            if not isinstance(calculator.wfs, FDWaveFunctions):
                raise RuntimeError(
                    'Linear response TDDFT supported only in real space mode')
            if calculator.wfs.kd.comm.size > 1:
                err_txt = 'Spin parallelization with Linear response '
                err_txt += "TDDFT. Use parallel = {'domain' : 'domain_only'} "
                err_txt += 'calculator parameter.'
                raise NotImplementedError(err_txt)
            if self.xc == 'GS':
                self.xc = calculator.hamiltonian.xc.name
            if calculator.input_parameters.mode != 'lcao':
                calculator.converge_wave_functions()
            if calculator.density.nct_G is None:
                spos_ac = calculator.initialize_positions()
                calculator.wfs.initialize(calculator.density,
                                          calculator.hamiltonian, spos_ac)

            self.update(calculator)

    def set(self, **kwargs):

        defaults = {
            'nspins': None,
            'eps': 0.001,
            'istart': 0,
            'jend': sys.maxsize,
            'energy_range': None,
            'xc': 'GS',
            'derivative_level': 1,
            'numscale': 0.00001,
            'txt': None,
            'filename': None,
            'finegrid': 2,
            'force_ApmB': False,  # for tests
            'eh_comm': None  # parallelization over eh-pairs
        }

        changed = False
        for key, value in defaults.items():
            if hasattr(self, key):
                value = getattr(self, key)  # do not overwrite
            setattr(self, key, kwargs.pop(key, value))
            if value != getattr(self, key):
                changed = True

        for key in kwargs:
            raise KeyError('Unknown key ' + key)

        return changed

    def set_calculator(self, calculator):
        self.calculator = calculator
#        self.force_ApmB = parameters['force_ApmB']
        self.force_ApmB = None  # XXX

    def analyse(self, what=None, out=None, min=0.1):
        """Print info about the transitions.

        Parameters:
          1. what: I list of excitation indicees, None means all
          2. out : I where to send the output, None means sys.stdout
          3. min : I minimal contribution to list (0<min<1)
        """
        if what is None:
            what = range(len(self))
        elif isinstance(what, int):
            what = [what]

        if out is None:
            out = sys.stdout

        for i in what:
            print(str(i) + ':', self[i].analyse(min=min), file=out)

    def update(self, calculator=None, **kwargs):

        changed = self.set(**kwargs)
        if calculator is not None:
            changed = True
            self.set_calculator(calculator)

        if not changed:
            return

        self.forced_update()

    def forced_update(self):
        """Recalc yourself."""
        if not self.force_ApmB:
            Om = OmegaMatrix
            name = 'LrTDDFT'
            if self.xc:
                xc = XC(self.xc)
                if hasattr(xc, 'hybrid') and xc.hybrid > 0.0:
                    Om = ApmB
                    name = 'LrTDDFThyb'
        else:
            Om = ApmB
            name = 'LrTDDFThyb'

        self.kss = KSSingles(calculator=self.calculator,
                             nspins=self.nspins,
                             eps=self.eps,
                             istart=self.istart,
                             jend=self.jend,
                             energy_range=self.energy_range,
                             txt=self.txt)

        self.Om = Om(self.calculator, self.kss,
                     self.xc, self.derivative_level, self.numscale,
                     finegrid=self.finegrid, eh_comm=self.eh_comm,
                     txt=self.txt)
        self.name = name

    def diagonalize(self, istart=None, jend=None,
                    energy_range=None, TDA=False):
        self.timer.start('diagonalize')
        self.timer.start('omega')
        self.Om.diagonalize(istart, jend, energy_range, TDA)
        self.timer.stop('omega')

        # remove old stuff
        self.timer.start('clean')
        while len(self):
            self.pop()
        self.timer.stop('clean')

        print('LrTDDFT digonalized:', file=self.txt)
        self.timer.start('build')
        for j in range(len(self.Om.kss)):
            self.append(LrTDDFTExcitation(self.Om, j))
            print(' ', str(self[-1]), file=self.txt)
        self.timer.stop('build')
        self.timer.stop('diagonalize')

    def get_Om(self):
        return self.Om

    def read(self, filename=None, fh=None):
        """Read myself from a file"""

        if fh is None:
            if filename.endswith('.gz'):
                try:
                    import gzip
                    f = gzip.open(filename)
                except:
                    f = open(filename, 'r')
            else:
                f = open(filename, 'r')
            self.filename = filename
        else:
            f = fh
            self.filename = None

        # get my name
        s = f.readline().replace('\n', '')
        self.name = s.split()[1]

        self.xc = f.readline().replace('\n', '').split()[0]
        values = f.readline().split()
        self.eps = float(values[0])
        if len(values) > 1:
            self.derivative_level = int(values[1])
            self.numscale = float(values[2])
            self.finegrid = int(values[3])
        else:
            # old writing style, use old defaults
            self.numscale = 0.001

        self.kss = KSSingles(filehandle=f)
        if self.name == 'LrTDDFT':
            self.Om = OmegaMatrix(kss=self.kss, filehandle=f,
                                  txt=self.txt)
        else:
            self.Om = ApmB(kss=self.kss, filehandle=f,
                           txt=self.txt)
        self.Om.Kss(self.kss)

        # check if already diagonalized
        p = f.tell()
        s = f.readline()
        if s != '# Eigenvalues\n':
            # go back to previous position
            f.seek(p)
        else:
            # load the eigenvalues
            n = int(f.readline().split()[0])
            for i in range(n):
                self.append(LrTDDFTExcitation(string=f.readline()))
            # load the eigenvectors
            f.readline()
            for i in range(n):
                values = f.readline().split()
                weights = [float(val) for val in values]
                self[i].f = np.array(weights)
                self[i].kss = self.kss

        if fh is None:
            f.close()

        # update own variables
        self.istart = self.Om.fullkss.istart
        self.jend = self.Om.fullkss.jend

    def singlets_triplets(self):
        """Split yourself into a singlet and triplet object"""

        slr = LrTDDFT(None, nspins=self.nspins, eps=self.eps,
                      istart=self.istart, jend=self.jend, xc=self.xc,
                      derivative_level=self.derivative_level,
                      numscale=self.numscale)
        tlr = LrTDDFT(None, nspins=self.nspins, eps=self.eps,
                      istart=self.istart, jend=self.jend, xc=self.xc,
                      derivative_level=self.derivative_level,
                      numscale=self.numscale)
        slr.Om, tlr.Om = self.Om.singlets_triplets()
        for lr in [slr, tlr]:
            lr.kss = lr.Om.fullkss
        return slr, tlr

    def single_pole_approximation(self, i, j):
        """Return the excitation according to the
        single pole approximation. See e.g.:
        Grabo et al, Theochem 501 (2000) 353-367
        """
        for ij, kss in enumerate(self.kss):
            if kss.i == i and kss.j == j:
                return sqrt(self.Om.full[ij][ij]) * Hartree
                return self.Om.full[ij][ij] / kss.energy * Hartree

    def __str__(self):
        string = ExcitationList.__str__(self)
        string += '# derived from:\n'
        string += self.Om.kss.__str__()
        return string

    def write(self, filename=None, fh=None):
        """Write current state to a file.

        'filename' is the filename. If the filename ends in .gz,
        the file is automatically saved in compressed gzip format.

        'fh' is a filehandle. This can be used to write into already
        opened files.
        """
        if mpi.rank == mpi.MASTER:
            if fh is None:
                if filename.endswith('.gz'):
                    try:
                        import gzip
                        f = gzip.open(filename, 'wb')
                    except:
                        f = open(filename, 'w')
                else:
                    f = open(filename, 'w')
            else:
                f = fh

            f.write('# ' + self.name + '\n')
            xc = self.xc
            if xc is None:
                xc = 'RPA'
            if self.calculator is not None:
                xc += ' ' + self.calculator.get_xc_functional()
            f.write(xc + '\n')
            f.write('%g %d %g %d' % (self.eps, int(self.derivative_level),
                                     self.numscale, int(self.finegrid)) + '\n')
            self.kss.write(fh=f)
            self.Om.write(fh=f)

            if len(self):
                f.write('# Eigenvalues\n')
                istart = self.istart
                if istart is None:
                    istart = self.kss.istart
                jend = self.jend
                if jend is None:
                    jend = self.kss.jend
                f.write('%d %d %d' % (len(self), istart, jend) + '\n')
                for ex in self:
                    f.write(ex.outstring())
                f.write('# Eigenvectors\n')
                for ex in self:
                    for w in ex.f:
                        f.write('%g ' % w)
                    f.write('\n')

            if fh is None:
                f.close()


def d2Excdnsdnt(dup, ddn):
    """Second derivative of Exc polarised"""
    res = [[0, 0], [0, 0]]
    for ispin in range(2):
        for jspin in range(2):
            res[ispin][jspin] = np.zeros(dup.shape)
            _gpaw.d2Excdnsdnt(dup, ddn, ispin, jspin, res[ispin][jspin])
    return res


def d2Excdn2(den):
    """Second derivative of Exc unpolarised"""
    res = np.zeros(den.shape)
    _gpaw.d2Excdn2(den, res)
    return res


class LrTDDFTExcitation(Excitation):

    def __init__(self, Om=None, i=None,
                 e=None, m=None, string=None):

        if string is not None:
            self.fromstring(string)
            return None

        # multiplicity comes from Kohn-Sham contributions
        self.fij = 1

        # define from the diagonalized Omega matrix
        if Om is not None:
            if i is None:
                raise RuntimeError

            ev = Om.eigenvalues[i]
            if ev < 0:
                # we reached an instability, mark it with a negative value
                self.energy = -sqrt(-ev)
            else:
                self.energy = sqrt(ev)
            self.f = Om.eigenvectors[i]
            self.kss = Om.kss

            self.kss.set_arrays()
            self.me = np.dot(self.f, self.kss.me)
            erat_k = np.sqrt(self.kss.energies / self.energy)
            wght_k = np.sqrt(self.kss.fij) * self.f
            ew_k = erat_k * wght_k
            self.mur = np.dot(ew_k, self.kss.mur)
            if self.kss.muv is not None:
                self.muv = np.dot(ew_k, self.kss.muv)
            else:
                self.muv = None
            if self.kss.magn is not None:
                self.magn = np.dot(1. / ew_k, self.kss.magn)
            else:
                self.magn = None

            return

        # define from energy and matrix element
        if e is not None:
            self.energy = e
            if m is None:
                if mur is None or muv is None or magn is None:
                    raise RuntimeError
                self.mur = mur
                self.muv = muv
                self.magn = magn
            else:
                self.me = m
            return

        raise RuntimeError

    def density_change(self, paw):
        """get the density change associated with this transition"""
        raise NotImplementedError

    def fromstring(self, string):
        l = string.split()
        self.energy = float(l.pop(0))
        if len(l) == 3:  # old writing style
            self.me = np.array([float(l.pop(0)) for i in range(3)])
        else:
            self.mur = np.array([float(l.pop(0)) for i in range(3)])
            self.me = - self.mur * sqrt(self.energy)
            self.muv = np.array([float(l.pop(0)) for i in range(3)])
            self.magn = np.array([float(l.pop(0)) for i in range(3)])

    def outstring(self):
        str = '%g ' % self.energy
        str += '  '
        for m in self.mur:
            str += '%12.4e' % m
        str += '  '
        for m in self.muv:
            str += '%12.4e' % m
        str += '  '
        for m in self.magn:
            str += '%12.4e' % m
        str += '\n'
        return str

    def __str__(self):
        m2 = np.sum(self.me * self.me)
        m = sqrt(m2)
        if m > 0:
            me = self.me / m
        else:
            me = self.me
        str = '<LrTDDFTExcitation> om=%g[eV] |me|=%g (%.2f,%.2f,%.2f)' % \
              (self.energy * Hartree, m, me[0], me[1], me[2])
        return str

    def analyse(self, min=.1):
        """Return an analysis string of the excitation"""
        osc = self.get_oscillator_strength()
        s = ('E=%.3f' % (self.energy * Hartree) + ' eV, ' +
             'f=%.5g' % osc[0] + ', (%.5g,%.5g,%.5g) ' %
             (osc[1], osc[2], osc[3]) + '\n')
           #'R=%.5g' % self.get_rotatory_strength() + ' cgs\n')

        def sqr(x):
            return x * x
        spin = ['u', 'd']
        min2 = sqr(min)
        rest = np.sum(self.f ** 2)
        for f, k in zip(self.f, self.kss):
            f2 = sqr(f)
            if f2 > min2:
                s += '  %d->%d ' % (k.i, k.j) + spin[k.pspin] + ' '
                s += '%.3g \n' % f2
                rest -= f2
        s += '  rest=%.3g' % rest
        return s


def photoabsorption_spectrum(excitation_list, spectrum_file=None,
                             e_min=None, e_max=None, delta_e=None,
                             folding='Gauss', width=0.1, comment=None):
    """Uniform absorption spectrum interface

    Parameters:
    ================= ===================================================
    ``exlist``        ExcitationList
    ``spectrum_file`` File name for the output file, STDOUT if not given
    ``e_min``         min. energy, set to cover all energies if not given
    ``e_max``         max. energy, set to cover all energies if not given
    ``delta_e``       energy spacing
    ``energyunit``    Energy unit, default 'eV'
    ``folding``       Gauss (default) or Lorentz
    ``width``         folding width in terms of the chosen energyunit
    ================= ===================================================
    all energies in [eV]
    """

    spectrum(exlist=excitation_list, filename=spectrum_file,
             emin=e_min, emax=e_max,
             de=delta_e, energyunit='eV',
             folding=folding, width=width,
             comment=comment)
