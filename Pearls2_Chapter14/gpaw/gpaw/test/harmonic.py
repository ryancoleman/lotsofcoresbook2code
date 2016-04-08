from math import pi

import numpy as np
from ase import Atoms
from ase.units import Hartree, Bohr

from gpaw import GPAW
from gpaw.xc import XC
from gpaw.test import equal
from gpaw.xc.kernel import XCNull
from gpaw.poisson import NoInteractionPoissonSolver


a = 4.0
x = Atoms(cell=(a, a, a))  # no atoms

class HarmonicPotential:
    def __init__(self, alpha):
        self.alpha = alpha
        self.vext_g = None

    def get_potential(self, gd):
        if self.vext_g is None:
            r_vg = gd.get_grid_point_coordinates()
            self.vext_g = self.alpha * ((r_vg - a / Bohr / 2)**2).sum(0)
        return self.vext_g


calc = GPAW(charge=-8,
            nbands=4,
            h=0.2,
            xc=XC(XCNull()),
            external=HarmonicPotential(0.5),
            poissonsolver=NoInteractionPoissonSolver(),
            eigensolver='cg')

x.calc = calc
x.get_potential_energy()

eigs = calc.get_eigenvalues()
equal(eigs[0], 1.5 * Hartree, 0.002)
equal(abs(eigs[1:] - 2.5 * Hartree).max(), 0, 0.003)
