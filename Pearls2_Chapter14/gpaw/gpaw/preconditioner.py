# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

from math import pi

import numpy as np

from gpaw.transformers import Transformer
from gpaw.fd_operators import Laplace

from gpaw.utilities.blas import axpy


class Preconditioner:
    def __init__(self, gd0, kin0, dtype=float, block=1):
        gd1 = gd0.coarsen()
        gd2 = gd1.coarsen()
        self.kin0 = kin0
        self.kin1 = Laplace(gd1, -0.5, 1, dtype)
        self.kin2 = Laplace(gd2, -0.5, 1, dtype)
        self.scratch0 = gd0.zeros((2, block), dtype, False)
        self.scratch1 = gd1.zeros((3, block),dtype, False)
        self.scratch2 = gd2.zeros((3, block), dtype, False)
        self.step = 0.66666666 / kin0.get_diagonal_element()

        self.restrictor_object0 = Transformer(gd0, gd1, 1,dtype)
        self.restrictor_object1 = Transformer(gd1, gd2, 1, dtype)
        self.interpolator_object2 = Transformer(gd2, gd1, 1, dtype)
        self.interpolator_object1 = Transformer(gd1, gd0, 1, dtype)
        self.restrictor0 = self.restrictor_object0.apply
        self.restrictor1 = self.restrictor_object1.apply
        self.interpolator2 = self.interpolator_object2.apply
        self.interpolator1 = self.interpolator_object1.apply

    def calculate_kinetic_energy(self, psit_xG, kpt):
        return None
        
    def __call__(self, residuals, kpt, ekin=None):
        nb = len(residuals) # number of bands
        phases = kpt.phase_cd
        step = self.step
        d0, q0 = self.scratch0[:,:nb]
        r1, d1, q1 = self.scratch1[:, :nb]
        r2, d2, q2 = self.scratch2[:, :nb]
        self.restrictor0(-residuals, r1, phases)
        d1[:] = 4 * step * r1
        self.kin1.apply(d1, q1, phases)
        q1 -= r1
        self.restrictor1(q1, r2, phases)
        d2 = 16 * step * r2
        self.kin2.apply(d2, q2, phases)
        q2 -= r2
        d2 -= 16 * step * q2
        self.interpolator2(d2, q1, phases)
        d1 -= q1
        self.kin1.apply(d1, q1, phases)
        q1 -= r1
        d1 -= 4 * step * q1
        self.interpolator1(-d1, d0, phases)
        self.kin0.apply(d0, q0, phases)
        q0 -= residuals
        axpy(-step, q0, d0)  # d0 -= step * q0
        d0 *= -1.0
        return d0

