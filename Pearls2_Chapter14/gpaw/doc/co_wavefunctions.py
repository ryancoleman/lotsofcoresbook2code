#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import Numeric as num
from gpaw import GPAW
from gpaw.spherical_harmonics import Y

a = 6.0
c = a / 2
d = 1.13
paw = GPAW('co.gpw', txt=None)

import pylab as p
dpi = 2*80
p.figure(figsize=(4, 3), dpi=dpi)
p.axes([0.15, 0.15, 0.8, 0.8])
import sys
if len(sys.argv) == 1:
    N = 1
else:
    N = int(sys.argv[1])
psit = paw.kpt_u[0].psit_nG[N]
psit = psit[:, 0, 0]
ng = len(psit)
x = num.arange(ng) * a / ng - c
p.plot(x, psit, 'bx', mew=2, label=r'$\tilde{\psi}$')

C = 'g'
for n in paw.nuclei:
    s = n.setup
    phi_j, phit_j = s.get_partial_waves()[:2]
    print(s.rcut_j[0])
    r = num.arange(30) * s.rcut_j[0] / 30
    phi_i = num.empty((s.ni, 59), num.Float)
    phit_i = num.empty((s.ni, 59), num.Float)
    x = num.empty(59, num.Float)
    x[29:] = r
    x[29::-1] = -r
    x *= paw.a0
    i = 0
    nj = len(phi_j)
    for j in range(nj):
        f = phi_j[j]
        ft = phit_j[j]
        l = f.get_angular_momentum_number()
        f = num.array([f(R) for R in r]) * r**l
        ft = num.array([ft(R) for R in r]) * r**l
        for m in range(2 * l + 1):
            L = l**2 + m
            phi_i[i + m, 29:] = f * Y(L, 1, 0, 0)
            phi_i[i + m, 29::-1] = f * Y(L, -1, 0, 0)
            phit_i[i + m, 29:] = ft * Y(L, 1, 0, 0)
            phit_i[i + m, 29::-1] = ft * Y(L, -1, 0, 0)
        i += 2 * l + 1
    assert i == s.ni
    P_i = n.P_uni[0, N]
    X = n.spos_c[0] * a - c
    p.plot(X + x, num.dot(P_i, phit_i), C + '-', lw=1,
           label=r'$\tilde{\psi}^%s$' % s.symbol)
    p.plot(X + x, num.dot(P_i, phi_i), C + '-', lw=2,
           label=r'$\psi^%s$' % s.symbol)
    C = 'r'
p.plot([-d / 2], [0], 'go', ms=2*0.8*4*80/a*1.0*0.53, mfc=None, label='_nolegend_')
p.plot([d / 2], [0], 'ro', ms=2*0.8*4*80/a*1.2*0.53, mfc=None, label='_nolegend_')
p.legend(loc='best')
p.xlabel(u'x [Å]')
p.ylabel(r'$\psi$')
#p.show()
p.savefig('co_wavefunctions.png', dpi=dpi)
