from __future__ import print_function, division

import functools
import os
import sys
from math import pi
from time import ctime

import numpy as np
from ase.units import Hartree
from ase.utils import devnull
from ase.utils.timing import timer, Timer
from scipy.special.orthogonal import p_roots

import gpaw.mpi as mpi
from gpaw import GPAW
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.response.chi0 import Chi0
from gpaw.response.wstc import WignerSeitzTruncatedCoulomb
from gpaw.wavefunctions.pw import PWDescriptor, count_reciprocal_vectors


def rpa(filename, ecut=200.0, blocks=1, extrapolate=4):
    """Calculate RPA energy.
    
    filename: str
        Name of restart-file.
    ecut: float
        Plane-wave cutoff.
    blocks: int
        Split polarizability matrix in this many blocks.
    extrapolate: int
        Number of cutoff energies to use for extrapolation.
    """
    name, ext = filename.rsplit('.', 1)
    assert ext == 'gpw'
    from gpaw.xc.rpa import RPACorrelation
    rpa = RPACorrelation(name, name + '-rpa.dat',
                         nblocks=blocks,
                         wstc=True,
                         txt=name + '-rpa.txt')
    rpa.calculate(ecut=ecut * (1 + 0.5 * np.arange(extrapolate))**(-2 / 3))

    
class RPACorrelation:
    def __init__(self, calc, xc='RPA', filename=None,
                 skip_gamma=False, qsym=True, nlambda=None,
                 nfrequencies=16, frequency_max=800.0, frequency_scale=2.0,
                 frequencies=None, weights=None,
                 world=mpi.world, nblocks=1, wstc=False,
                 txt=sys.stdout):

        if isinstance(calc, str):
            calc = GPAW(calc, txt=None, communicator=mpi.serial_comm)
        self.calc = calc

        if world.rank != 0:
            txt = devnull
        elif isinstance(txt, str):
            txt = open(txt, 'w')
        self.fd = txt

        self.timer = Timer()
        
        if frequencies is None:
            frequencies, weights = get_gauss_legendre_points(nfrequencies,
                                                             frequency_max,
                                                             frequency_scale)
            user_spec = False
        else:
            assert weights is not None
            user_spec = True
            
        self.omega_w = frequencies / Hartree
        self.weight_w = weights / Hartree

        if nblocks > 1:
            assert len(self.omega_w) % nblocks == 0
            assert wstc

        self.wstc = wstc
        self.nblocks = nblocks
        self.world = world

        self.skip_gamma = skip_gamma
        self.ibzq_qc = None
        self.weight_q = None
        self.initialize_q_points(qsym)

        # Energies for all q-vetors and cutoff energies:
        self.energy_qi = []

        self.filename = filename

        self.print_initialization(xc, frequency_scale, nlambda, user_spec)

    def initialize_q_points(self, qsym):
        kd = self.calc.wfs.kd
        self.bzq_qc = kd.get_bz_q_points(first=True)

        if not qsym:
            self.ibzq_qc = self.bzq_qc
            self.weight_q = np.ones(len(self.bzq_qc)) / len(self.bzq_qc)
        else:
            U_scc = kd.symmetry.op_scc
            self.ibzq_qc = kd.get_ibz_q_points(self.bzq_qc, U_scc)[0]
            self.weight_q = kd.q_weights

    def read(self):
        lines = open(self.filename).readlines()[1:]
        n = 0
        self.energy_qi = []
        nq = len(lines) // len(self.ecut_i)
        for q_c in self.ibzq_qc[:nq]:
            self.energy_qi.append([])
            for ecut in self.ecut_i:
                q1, q2, q3, ec, energy = [float(x)
                                          for x in lines[n].split()]
                self.energy_qi[-1].append(energy / Hartree)
                n += 1

                if (abs(q_c - (q1, q2, q3)).max() > 1e-4 or
                    abs(int(ecut * Hartree) - ec) > 0):
                    self.energy_qi = []
                    return

        print('Read %d q-points from file: %s' % (nq, self.filename),
              file=self.fd)
        print(file=self.fd)

    def write(self):
        if self.world.rank == 0 and self.filename:
            fd = open(self.filename, 'w')
            print('#%9s %10s %10s %8s %12s' %
                  ('q1', 'q2', 'q3', 'E_cut', 'E_c(q)'), file=fd)
            for energy_i, q_c in zip(self.energy_qi, self.ibzq_qc):
                for energy, ecut in zip(energy_i, self.ecut_i):
                    print('%10.4f %10.4f %10.4f %8d   %r' %
                          (tuple(q_c) + (ecut * Hartree, energy * Hartree)),
                          file=fd)

    def calculate(self, ecut, nbands=None, spin=False):
        """Calculate RPA correlation energy for one or several cutoffs.

        ecut: float or list of floats
            Plane-wave cutoff(s).
        nbands: int
            Number of bands (defaults to number of plane-waves).
        spin: bool
            Separate spin in response funtion.
            (Only needed for beyond RPA methods that inherit this function).
        """

        p = functools.partial(print, file=self.fd)

        if isinstance(ecut, (float, int)):
            ecut = ecut * (1 + 0.5 * np.arange(6))**(-2 / 3)
        self.ecut_i = np.asarray(np.sort(ecut)) / Hartree
        ecutmax = max(self.ecut_i)

        if nbands is None:
            p('Response function bands : Equal to number of plane waves')
        else:
            p('Response function bands : %s' % nbands)
        p('Plane wave cutoffs (eV) :', end='')
        for e in self.ecut_i:
            p(' {0:.3f}'.format(e * Hartree), end='')
        p()
        p()

        if self.filename and os.path.isfile(self.filename):
            self.read()
            self.world.barrier()

        chi0 = Chi0(self.calc, 1j * Hartree * self.omega_w, eta=0.0,
                    intraband=False, hilbert=False,
                    txt=self.fd, timer=self.timer, world=self.world,
                    no_optical_limit=self.wstc,
                    nblocks=self.nblocks)

        self.blockcomm = chi0.blockcomm
        
        wfs = self.calc.wfs
        
        if self.wstc:
            with self.timer('WSTC-init'):
                p('Using Wigner-Seitz truncated Coulomb potential.')
                self.wstc = WignerSeitzTruncatedCoulomb(
                    wfs.gd.cell_cv, wfs.kd.N_c, self.fd)

        nq = len(self.energy_qi)
        nw = len(self.omega_w)
        nGmax = max(count_reciprocal_vectors(ecutmax, wfs.gd, q_c)
                    for q_c in self.ibzq_qc[nq:])
        mynGmax = (nGmax + self.nblocks - 1) // self.nblocks
        
        nx = (1 + spin) * nw * mynGmax * nGmax
        A1_x = np.empty(nx, complex)
        if self.nblocks > 1:
            A2_x = np.empty(nx, complex)
        else:
            A2_x = None
        
        self.timer.start('RPA')
        
        for q_c in self.ibzq_qc[nq:]:
            if np.allclose(q_c, 0.0) and self.skip_gamma:
                self.energy_qi.append(len(self.ecut_i) * [0.0])
                self.write()
                p('Not calculating E_c(q) at Gamma')
                p()
                continue

            thisqd = KPointDescriptor([q_c])
            pd = PWDescriptor(ecutmax, wfs.gd, complex, thisqd)
            nG = pd.ngmax
            mynG = (nG + self.nblocks - 1) // self.nblocks
            chi0.Ga = self.blockcomm.rank * mynG
            chi0.Gb = min(chi0.Ga + mynG, nG)
            
            shape = (1 + spin, nw, chi0.Gb - chi0.Ga, nG)
            chi0_swGG = A1_x[:np.prod(shape)].reshape(shape)
            chi0_swGG[:] = 0.0
            
            if self.wstc or np.allclose(q_c, 0.0):
                # Wings (x=0,1) and head (G=0) for optical limit and three
                # directions (v=0,1,2):
                chi0_swxvG = np.zeros((1 + spin, nw, 2, 3, nG), complex)
                chi0_swvv = np.zeros((1 + spin, nw, 3, 3), complex)
            else:
                chi0_swxvG = None
                chi0_swvv = None

            Q_aGii = chi0.initialize_paw_corrections(pd)

            # First not completely filled band:
            m1 = chi0.nocc1
            p('# %s  -  %s' % (len(self.energy_qi), ctime().split()[-2]))
            p('q = [%1.3f %1.3f %1.3f]' % tuple(q_c))

            energy_i = []
            for ecut in self.ecut_i:
                if ecut == ecutmax:
                    # Nothing to cut away:
                    cut_G = None
                    m2 = nbands or nG
                else:
                    cut_G = np.arange(nG)[pd.G2_qG[0] <= 2 * ecut]
                    m2 = len(cut_G)

                p('E_cut = %d eV / Bands = %d:' % (ecut * Hartree, m2))
                self.fd.flush()

                energy = self.calculate_q(chi0, pd,
                                          chi0_swGG, chi0_swxvG, chi0_swvv,
                                          Q_aGii, m1, m2, cut_G, A2_x)
                energy_i.append(energy)
                m1 = m2

                a = 1 / chi0.kncomm.size
                if ecut < ecutmax and a != 1.0:
                    # Chi0 will be summed again over chicomm, so we divide
                    # by its size:
                    chi0_swGG *= a
                    if chi0_swxvG is not None:
                        chi0_swxvG *= a
                        chi0_swvv *= a

            self.energy_qi.append(energy_i)
            self.write()
            p()

        e_i = np.dot(self.weight_q, np.array(self.energy_qi))
        p('==========================================================')
        p()
        p('Total correlation energy:')
        for e_cut, e in zip(self.ecut_i, e_i):
            p('%6.0f:   %6.4f eV' % (e_cut * Hartree, e * Hartree))
        p()

        self.energy_qi = []  # important if another calculation is performed

        if len(e_i) > 1:
            self.extrapolate(e_i)

        p('Calculation completed at: ', ctime())
        p()

        self.timer.stop('RPA')
        self.timer.write(self.fd)
        self.fd.flush()
        
        return e_i * Hartree

    @timer('chi0(q)')
    def calculate_q(self, chi0, pd,
                    chi0_swGG, chi0_swxvG, chi0_swvv, Q_aGii, m1, m2, cut_G,
                    A2_x):
        chi0_wGG = chi0_swGG[0]
        if chi0_swxvG is not None:
            chi0_wxvG = chi0_swxvG[0]
            chi0_wvv = chi0_swvv[0]
        else:
            chi0_wxvG = None
            chi0_wvv = None
        chi0._calculate(pd, chi0_wGG, chi0_wxvG, chi0_wvv,
                        Q_aGii, m1, m2, [0, 1])

        print('E_c(q) = ', end='', file=self.fd)

        chi0_wGG = chi0.redistribute(chi0_wGG, A2_x)
        
        if not pd.kd.gamma or self.wstc:
            e = self.calculate_energy(pd, chi0_wGG, cut_G)
            print('%.3f eV' % (e * Hartree), file=self.fd)
            self.fd.flush()
        else:
            e = 0.0
            for v in range(3):
                chi0_wGG[:, 0] = chi0_wxvG[:, 0, v]
                chi0_wGG[:, :, 0] = chi0_wxvG[:, 1, v]
                chi0_wGG[:, 0, 0] = chi0_wvv[:, v, v]
                ev = self.calculate_energy(pd, chi0_wGG, cut_G)
                e += ev
                print('%.3f' % (ev * Hartree), end='', file=self.fd)
                if v < 2:
                    print('/', end='', file=self.fd)
                else:
                    print(' eV', file=self.fd)
                    self.fd.flush()
            e /= 3

        return e

    @timer('Energy')
    def calculate_energy(self, pd, chi0_wGG, cut_G):
        """Evaluate correlation energy from chi0."""

        if self.wstc:
            invG_G = (self.wstc.get_potential(pd) / (4 * pi))**0.5
        else:
            G_G = pd.G2_qG[0]**0.5  # |G+q|
            if pd.kd.gamma:
                G_G[0] = 1.0
            invG_G = 1.0 / G_G
            
        if cut_G is not None:
            invG_G = invG_G[cut_G]

        nG = len(invG_G)

        e_w = []
        for chi0_GG in chi0_wGG:
            if cut_G is not None:
                chi0_GG = chi0_GG.take(cut_G, 0).take(cut_G, 1)

            e_GG = (np.eye(nG) -
                    4 * np.pi * chi0_GG * invG_G * invG_G[:, np.newaxis])
            e = np.log(np.linalg.det(e_GG)) + nG - np.trace(e_GG)
            e_w.append(e.real)

        E_w = np.zeros_like(self.omega_w)
        self.blockcomm.all_gather(np.array(e_w), E_w)
        energy = np.dot(E_w, self.weight_w) / (2 * np.pi)
        self.E_w = E_w
        return energy

    def extrapolate(self, e_i):
        print('Extrapolated energies:', file=self.fd)
        ex_i = []
        for i in range(len(e_i) - 1):
            e1, e2 = e_i[i:i + 2]
            x1, x2 = self.ecut_i[i:i + 2]**-1.5
            ex = (e1 * x2 - e2 * x1) / (x2 - x1)
            ex_i.append(ex)

            print('  %4.0f -%4.0f:  %5.3f eV' % (self.ecut_i[i] * Hartree,
                                                 self.ecut_i[i + 1] * Hartree,
                                                 ex * Hartree),
                  file=self.fd)
        print(file=self.fd)
        self.fd.flush()

        return e_i * Hartree

    def print_initialization(self, xc, frequency_scale, nlambda, user_spec):
        p = functools.partial(print, file=self.fd)
        p('----------------------------------------------------------')
        p('Non-self-consistent %s correlation energy' % xc)
        p('----------------------------------------------------------')
        p('Started at:  ', ctime())
        p()
        p('Atoms                          :',
          self.calc.atoms.get_chemical_formula(mode='hill'))
        p('Ground state XC functional     :', self.calc.hamiltonian.xc.name)
        p('Valence electrons              :', self.calc.wfs.setups.nvalence)
        p('Number of bands                :', self.calc.wfs.bd.nbands)
        p('Number of spins                :', self.calc.wfs.nspins)
        p('Number of k-points             :', len(self.calc.wfs.kd.bzk_kc))
        p('Number of irreducible k-points :', len(self.calc.wfs.kd.ibzk_kc))
        p('Number of q-points             :', len(self.bzq_qc))
        p('Number of irreducible q-points :', len(self.ibzq_qc))
        p()
        for q, weight in zip(self.ibzq_qc, self.weight_q):
            p('    q: [%1.4f %1.4f %1.4f] - weight: %1.3f' %
              (q[0], q[1], q[2], weight))
        p()
        p('----------------------------------------------------------')
        p('----------------------------------------------------------')
        p()
        if nlambda is None:
            p('Analytical coupling constant integration')
        else:
            p('Numerical coupling constant integration using', nlambda,
              'Gauss-Legendre points')
        p()
        p('Frequencies')
        if not user_spec:
            p('    Gauss-Legendre integration with %s frequency points' %
              len(self.omega_w))
            p('    Transformed from [0,oo] to [0,1] using e^[-aw^(1/B)]')
            p('    Highest frequency point at %5.1f eV and B=%1.1f' %
              (self.omega_w[-1] * Hartree, frequency_scale))
        else:
            p('    User specified frequency integration with',
              len(self.omega_w), 'frequency points')
        p()
        p('Parallelization')
        p('    Total number of CPUs          : % s' % self.world.size)
        p('    G-vector decomposition       : % s' % self.nblocks)
        p('    K-point/band decomposition    : % s' %
          (self.world.size // self.nblocks))
        p()


def get_gauss_legendre_points(nw=16, frequency_max=800.0, frequency_scale=2.0):
    y_w, weights_w = p_roots(nw)
    y_w = y_w.real
    ys = 0.5 - 0.5 * y_w
    ys = ys[::-1]
    w = (-np.log(1 - ys))**frequency_scale
    w *= frequency_max / w[-1]
    alpha = (-np.log(1 - ys[-1]))**frequency_scale / frequency_max
    transform = (-np.log(1 - ys))**(frequency_scale - 1) \
        / (1 - ys) * frequency_scale / alpha
    return w, weights_w * transform / 2

    
description = 'Run RPA-correlation calculation.'


def main():
    import optparse
    parser = optparse.OptionParser(usage='Usage: %prog <gpw-file> [options]',
                                   description=description)
    add = parser.add_option
    
    add('-e', '--cut-off', type=float, default=100, meta='ECUT',
        help='Plane-wave cut off energy (eV) for polarization function.')
    add('-b', '--blocks', type=int, default=1, meta='N',
        help='Split polarization matrix in N blocks.')
    
    opts, args = parser.parse_args()
    if len(args) == 0:
        parser.error('No gpw-file!')
    if len(args) > 1:
        parser.error('Too many arguments!')
    
    name = args[0]
    assert name.endswith('.gpw')
    
    rpa = RPACorrelation(name,
                         txt=name[:-3] + 'rpa.txt',
                         wstc=True, nblocks=opts.blocks)
    rpa.calculate([opts.cut_off])


if __name__ == '__main__':
    main()
