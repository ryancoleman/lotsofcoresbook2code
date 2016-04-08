from __future__ import print_function
from math import pi, sin, cos
import numpy as np

from ase.units import Bohr
import _gpaw
from gpaw.poisson import PoissonSolver
from gpaw.pes.state import State
from gpaw.analyse.expandyl import AngularIntegral


class PlaneWave(State):

    """Plane wave state."""

    def get_grid(self, k_c, r0=None):
        if r0 is None:
            r0_c = np.array([0., 0., 0.])
        else:
            r0_c = np.array(r0)
        pw_G = self.gd.zeros(dtype=complex)
        python = True
        python = False
        if python:
            r_c = np.array([0., 0., 0.])
            gd = self.gd
            for i in range(gd.n_c[0]):
                r_c[0] = (gd.beg_c[0] + i) * gd.h_c[0]
                for j in range(gd.n_c[1]):
                    r_c[1] = (gd.beg_c[1] + j) * gd.h_c[1]
                    for k in range(gd.n_c[2]):
                        r_c[2] = (gd.beg_c[2] + k) * gd.h_c[2]
                        phase = np.dot(k_c, r_c - r0_c)
                        pw_G[i, j, k] = cos(phase) + 1j * sin(phase)
        else:
            gd = self.gd
            _gpaw.plane_wave_grid(gd.beg_c, gd.end_c,
                                  gd.h_cv.diagonal().copy(),
                                  k_c, r0_c, pw_G)
        pw_G /= (2 * pi) ** (3. / 2.)
        return pw_G


class ZerothOrder1(State):

    def __init__(self, calculator):
        self.calculator = calculator
        State.__init__(self, calculator.gd)
        self.pw = PlaneWave(self.gd)

        # get effective potential
        hamiltonian = self.calculator.hamiltonian
        if 1:
            self.vt_G = hamiltonian.vt_sG[0]  # XXX treat spin
        else:
            self.vt_G = np.where(hamiltonian.vt_sG[0] > 0,
                                 0.0, hamiltonian.vt_sG[0])
        self.intvt = self.gd.integrate(self.vt_G)
        print('# int(vt_G)=', self.intvt, np.sometrue(self.vt_G > 0))

        self.solve()

    def get_grid(self, k_c, r0):
        print('Correction:', self.gd.integrate(self.corrt_G))
        if not hasattr(self, 'written'):
            print('r0', r0, r0 * Bohr)
            dR = 0.3
            AI = AngularIntegral(r0 * Bohr, self.gd, dR=dR)
            r_R = AI.radii()
            psi_R = AI.average(self.corrt_G)
            v_R = AI.average(self.vt_G)
            f = open('radial_dR' + str(dR) + '.dat', 'w')
            print('# R  v(R)    psi(R)', file=f)
            for r, psi, v in zip(r_R, psi_R, v_R):
                print(r, psi, v, file=f)
            f.close()
            self.written = True

        return self.pw.get_grid(k_c, r0) - 1e6 * self.corrt_G

    def solve(self):
        hamiltonian = self.calculator.hamiltonian
        self.poisson = PoissonSolver(nn=hamiltonian.poisson.nn)
        self.poisson.set_grid_descriptor(self.gd)
        self.poisson.initialize()

        corrt_G = self.gd.empty()
        self.poisson.solve(corrt_G, self.vt_G, charge=None)
        corrt_G /= (2 * pi) ** (5. / 2)

        self.corrt_G = corrt_G
