# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""Grid transformers.

This module defines tools for doing interpolations/restrictions between
differentt uniform 3D grids.
"""

from __future__ import division
from __future__ import print_function
from math import pi
import numpy as np

from gpaw import debug
from gpaw.utilities import is_contiguous
import _gpaw


class _Transformer:
    def __init__(self, gdin, gdout, nn=1, dtype=float):
        self.gdin = gdin
        self.gdout = gdout
        self.nn = nn
        assert 1 <= nn <= 4
        self.dtype = dtype

        pad_cd = np.empty((3, 2), int)
        neighborpad_cd = np.empty((3, 2), int)
        skip_cd = np.empty((3, 2), int)

        if (gdin.N_c == 2 * gdout.N_c).all():
            # Restriction:
            pad_cd[:, 0] = 2 * nn - 1 - 2 * gdout.beg_c + gdin.beg_c
            pad_cd[:, 1] = 2 * nn - 2 + 2 * gdout.end_c - gdin.end_c
            neighborpad_cd[:, 0] = 2 * nn - 2 + 2 * gdout.beg_c - gdin.beg_c
            neighborpad_cd[:, 1] = 2 * nn - 1 - 2 * gdout.end_c + gdin.end_c
            self.interpolate = False
        else:
            assert (gdout.N_c == 2 * gdin.N_c).all()
            # Interpolation:
            pad_cd[:, 0] = nn - 1 - gdout.beg_c // 2 + gdin.beg_c
            pad_cd[:, 1] = nn + gdout.end_c // 2 - gdin.end_c
            neighborpad_cd[:, 0] = nn + gdout.beg_c // 2 - gdin.beg_c
            neighborpad_cd[:, 1] = nn - 1 - gdout.end_c // 2 + gdin.end_c
            skip_cd[:, 0] = gdout.beg_c % 2
            skip_cd[:, 1] = gdout.end_c % 2
            self.interpolate = True

            inpoints = (gdin.n_c[0] + 2 * nn - 1) * (gdin.n_c[1] + 2 * nn - 1)
            outpoints = gdout.n_c[0] * gdout.n_c[1]
            
            if inpoints > outpoints:
                points = ' x '.join([str(N) for N in gdin.N_c])
                raise ValueError('Cannot construct interpolator.  Grid %s '
                                 'may be too small' % points)

        assert np.alltrue(pad_cd.ravel() >= 0)
        self.ngpin = tuple(gdin.n_c)
        self.ngpout = tuple(gdout.n_c)
        assert dtype in [float, complex]

        self.pad_cd = pad_cd
        self.neighborpad_cd = neighborpad_cd
        self.skip_cd = skip_cd

        if gdin.comm.size > 1:
            comm = gdin.comm.get_c_object()
        else:
            comm = None
        
        self.transformer = _gpaw.Transformer(gdin.n_c, gdout.n_c, 
                                             2 * nn, pad_cd, 
                                             neighborpad_cd, skip_cd,
                                             gdin.neighbor_cd,
                                             dtype == float, comm,
                                             self.interpolate)
        
    def apply(self, input, output=None, phases=None):
        if output is None:
            output = self.gdout.empty(input.shape[:-3], dtype=self.dtype)
        self.transformer.apply(input, output, phases)
        return output

    def get_async_sizes(self):
        return self.transformer.get_async_sizes()


class TransformerWrapper:
    def __init__(self, transformer):
        self.transformer = transformer
        self.dtype = transformer.dtype
        self.ngpin = transformer.ngpin
        self.ngpout = transformer.ngpout
        self.nn = transformer.nn

    def apply(self, input, output=None, phases=None):
        assert is_contiguous(input, self.dtype)
        assert input.shape[-3:] == self.ngpin
        if output is not None:
            assert is_contiguous(output, self.dtype)
            assert output.shape[-3:] == self.ngpout
        assert (self.dtype == float or
                (phases.dtype == complex and
                 phases.shape == (3, 2)))

        return self.transformer.apply(input, output, phases)
        
    def get_async_sizes(self):
        return self.transformer.get_async_sizes()


def Transformer(gdin, gdout, nn=1, dtype=float):
    if nn != 9:
        t = _Transformer(gdin, gdout, nn, dtype)
        if debug:
            t = TransformerWrapper(t)
        return t
    class T:
        nn = 1
        def apply(self, input, output, phases=None):
            output[:] = input
    return T()


def multiple_transform_apply(transformerlist, inputs, outputs, phases=None):
    return _gpaw.multiple_transform_apply(transformerlist, inputs, outputs, 
                                          phases)


def coefs(k, p):
    for i in range(0, k * p, p):
        print('%2d' % i, end=' ')
        for x in range((k // 2 - 1) * p, k // 2 * p + 1):
            n = 1
            d = 1
            for j in range(0, k * p, p):
                if j == i:
                    continue
                n *= x - j
                d *= i - j
            print('%14.16f' % (n / d), end=' ')
        print()
