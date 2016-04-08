import os
import sys
import pickle

import matplotlib
#matplotlib.use('Agg')

import numpy as np
from ase.atoms import Atoms
from ase.utils import devnull
from ase.data import atomic_numbers, atomic_names, covalent_radii

from gpaw.atom.generator import Generator
from gpaw.atom.configurations import parameters
from gpaw.atom.analyse_setup import analyse
from gpaw import GPAW, ConvergenceError, Mixer, FermiDirac
import gpaw.mpi as mpi

b0 = {
      'H': 0.753,
      'He': 2.7,
      'Li': 2.99,
      'Be': 2.618,
      'B': 1.694,
      'C': 1.279,
      'N': 1.102,
      'O': 1.234,
      'F': 1.413,
      'Ne': 2.9,
      'Na': 3.289,
      'Mg': 3.5,
      'Al': 2.868,
      'Si': 2.218,
      'P': 1.878,
      'S': 1.893,
      'Cl': 1.989,
      'Ar': 3.7,
      'K': 4.108,
      'Ca': 2.805,
      'Sc': 2.3,
      'Ti': 2.0,
      'V': 1.82,
      'Cr': 1.75,
      'Mn': 1.78,
      'Fe': 1.850,
      'Co': 2.0,
      'Ni': 2.1,
      'Cu': 2.281,
      'Zn': 3.4,
      'Ga': 2.837,
      'Ge': 2.319,
      'As': 2.071,
      'Se': 2.154,
      'Br': 2.281,
      'Kr': 4.8,
      'Rb': 4.360,
      'Sr': 4.5,
      'Y': 2.67,
      'Zr': 2.36,
      'Nb': 2.14,
      'Mo': 1.95,
      'Ru': 2.1,
      'Rh': 2.20,
      'Pd': 2.485,
      'Ag': 2.626,
      'Cd': 3.6,
      'In': 3.1,
      'Sb': 2.5,
      'Xe': 4.5,
      'Cs': 4.819,
      'Ba': 4.60,
      'La': 2.872,
      'Hf': 2.38,
      'Ta': 2.2,
      'W': 2.1,
      'Re': 2.1,
      'Os': 2.2,
      'Ir': 2.227,
      'Pt': 2.373,
      'Au': 2.555,
      'Hg': 3.6,
      'Pb': 2.85,
      'Bi': 2.6,
      'Rn': 4.7,
      }

class MakeSetupPageData:
    def __init__(self, symbol):
        self.symbol = symbol
        self.Z = atomic_numbers[symbol]
        self.name = atomic_names[self.Z]
        self.parameters = dict(occupations=FermiDirac(width=0.1), xc='PBE')
        if mpi.rank == 0:
            self.log = sys.stdout
        else:
            self.log = devnull

    def run(self):
        if os.path.isfile(self.symbol + '.pckl'):
            self.log.write('Skipping %s\n' % self.symbol)
            return
        mpi.world.barrier()
        if mpi.rank == 0:
            self.file = open(self.symbol + '.pckl', 'w')

        self.generate_setup()
        self.prepare_box()
        self.eggbox()
        self.dimer()
        self.pickle()

    def generate_setup(self):
        if mpi.rank == 0:
            gen = Generator(self.symbol, 'PBE', scalarrel=True)
            gen.run(logderiv=True, **parameters[self.symbol])
            data = analyse(gen, show=False)

            g = np.arange(gen.N)
            r_g = gen.r
            dr_g = gen.beta * gen.N / (gen.N - g)**2
            rcutcomp = gen.rcutcomp
            rcutfilter = gen.rcutfilter

            # Find cutoff for core density:
            if gen.Nc == 0:
                rcore = 0.5
            else:
                N = 0.0
                g = gen.N - 1
                while N < 1e-7:
                    N += np.sqrt(4 * np.pi) * gen.nc[g] * r_g[g]**2 * dr_g[g]
                    g -= 1
                rcore = r_g[g]

            nlfer = []
            for j in range(gen.njcore):
                nlfer.append((gen.n_j[j], gen.l_j[j], gen.f_j[j], gen.e_j[j],
                              0.0))
            for n, l, f, eps in zip(gen.vn_j, gen.vl_j, gen.vf_j, gen.ve_j):
                nlfer.append((n, l, f, eps, gen.rcut_l[l]))

            self.data = dict(Z=gen.Z,
                             Nv=gen.Nv,
                             Nc=gen.Nc,
                             rcutcomp=rcutcomp,
                             rcutfilter=rcutfilter,
                             rcore=rcore,
                             Ekin=gen.Ekin,
                             Epot=gen.Epot,
                             Exc=gen.Exc,
                             nlfer=nlfer)

    def prepare_box(self):
        if symbol in b0:
            self.d0 = b0[symbol]
        else:
            self.d0 = 2 * covalent_radii[self.Z]

        if symbol in ['He', 'Ne', 'Ar', 'Kr']:
            self.a = round(2 / np.sqrt(3) * self.d0 / 0.2 / 4) * 4 * 0.2
        else:
            self.a = round(max(2.5 * self.d0, 5.5) / 0.2 / 4) * 4 * 0.2

        gmin = 4 * int(self.a / 0.22 / 4 + 0.5)
        gmax = 4 * int(self.a / 0.12 / 4 + 0.5)
        self.ng = (gmax + 4 - gmin) // 4
        self.gridspacings = self.a / np.arange(gmin, gmax + 4, 4)

    def eggbox(self):
        atom = Atoms(self.symbol, pbc=True, cell=(self.a, self.a, self.a))

        negg = 25
        self.Eegg = np.zeros((self.ng, negg))
        self.Fegg = np.zeros((self.ng, negg))

        eigensolver = 'rmm-diis'
        if self.symbol in ['Na', 'Mg']: eigensolver = 'cg'
        if self.symbol in ['Sc', 'Ti', 'V', 'Mn', 'Ni',
                           'Zn']: eigensolver = 'cg'
        if self.symbol in ['Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd',
                           'Ag']: eigensolver = 'cg'
        if self.symbol in ['Sn', 'Sb', 'Te', 'Ba']: eigensolver = 'cg'
        if self.symbol in ['Hf', 'Ta', 'W', 'Re', 'Os',
                           'Ir', 'Pt', 'Hg', 'Pb', 'Bi']: eigensolver = 'cg'
        for i in range(self.ng):
            h = self.gridspacings[i]
            calc = GPAW(h=h, txt='%s-eggbox-%.3f.txt' % (self.symbol, h),
                        mixer=Mixer(beta=0.1, nmaxold=5, weight=50),
                        eigensolver=eigensolver,
                        maxiter=300,
                        nbands=-10,
                        **self.parameters)
            atom.set_calculator(calc)

            for j in range(negg):
                x = h * j / (2 * negg - 2)
                atom[0].x = x
                try:
                    e = calc.get_potential_energy(atom, force_consistent=True)
                    self.Eegg[i, j] = e
                except ConvergenceError:
                    raise
                self.Fegg[i, j] = atom.get_forces()[0, 0]

    def dimer(self):
        dimer = Atoms([self.symbol, self.symbol],
                      pbc=True, cell=(self.a, self.a, self.a))

        self.Edimer = np.zeros((self.ng, 7))
        self.Fdimer = np.zeros((self.ng, 7, 2))

        q0 = self.d0 / np.sqrt(3)
        eigensolver = 'rmm-diis'
        if self.symbol in ['Na', 'Mg']: eigensolver = 'cg'
        if self.symbol in ['Sc', 'Ti', 'V', 'Mn', 'Ni',
                           'Zn']: eigensolver = 'cg'
        if self.symbol in ['Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd',
                           'Ag']: eigensolver = 'cg'
        if self.symbol in ['Sn', 'Sb', 'Te', 'Ba']: eigensolver = 'cg'
        if self.symbol in ['Hf', 'Ta', 'W', 'Re', 'Os',
                           'Ir', 'Pt', 'Hg', 'Pb', 'Bi']: eigensolver = 'cg'
        for i in range(self.ng):
            h = self.gridspacings[i]
            calc = GPAW(h=h, txt='%s-dimer-%.3f.txt' % (self.symbol, h),
                        mixer=Mixer(beta=0.1, nmaxold=5, weight=50),
                        #mixer=Mixer(beta=0.05, nmaxold=7, weight=100),
                        eigensolver=eigensolver,
                        maxiter=300,
                        nbands=-10,
                        **self.parameters)
            dimer.set_calculator(calc)

            y = []
            for j in range(-3, 4):
                q = q0 * (1 + j * 0.02)
                dimer.positions[1] = (q, q, q)
                try:
                    e = calc.get_potential_energy(dimer, force_consistent=True)
                    self.Edimer[i, j + 3] = e
                except ConvergenceError:
                    raise
                self.Fdimer[i, j + 3] = dimer.get_forces()[:, 0]

    def pickle(self):
        if mpi.rank == 0:
            self.data.update({'d0': self.d0,
                              'a': self.a,
                              'gridspacings': self.gridspacings,
                              'Eegg': self.Eegg,
                              'Fegg': self.Fegg,
                              'Edimer': self.Edimer,
                              'Fdimer': self.Fdimer})
            pickle.dump(self.data, self.file, pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) == 0:
        args = parameters.keys()
    for symbol in args:
        MakeSetupPageData(symbol).run()
