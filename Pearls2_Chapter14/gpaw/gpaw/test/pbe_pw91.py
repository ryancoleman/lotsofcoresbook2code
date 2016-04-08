from math import pi
import numpy as np
from gpaw.xc.libxc import LibXC

x0 = LibXC('LDA_X')

def f0(xc, rs, s):
    n = 3 / (4 * pi * rs**3)
    third = 1.0 / 3.0
    kF = (3 * pi**2 * n)**third
    a2 = (2 * kF * n * s)**2
    e = np.zeros(1)
    xc.calculate(e,
                 np.array([[n]]), np.zeros((1, 1)),
                 np.array([[a2]]), np.zeros((1, 1)))
    exc = n * e[0]
    x0.calculate(e, np.array([[n]]), np.zeros((1, 1)))
    ex0 = n * e[0]
    return exc / ex0

def f1(xc, rs, s):
    n = 3 / (4 * pi * rs**3)
    na = 2 * n
    third = 1.0 / 3.0
    kF = (3 * pi**2 * n)**third
    rsa = (3 / pi / 4 / na)**third
    a2 = (2 * kF * n * s)**2
    e = np.zeros(1)
    xc.calculate(e,
                 np.array([[n], [0]]), np.zeros((2, 1)),
                 np.array([[a2], [0], [0]]), np.zeros((3, 1)))
    exc = n * e[0]
    x0.calculate(e, np.array([[n]]), np.zeros((1, 1)))
    ex0 = n * e[0]
    return exc / ex0

pbe = LibXC('PBE')
pw91 = LibXC('PW91')
assert abs(f0(pbe, 2, 3) - 1.58) < 0.01
assert abs(f1(pbe, 2, 3) - 1.88) < 0.01
assert abs(f0(pw91, 2, 3) - 1.60) < 0.01
assert abs(f1(pw91, 2, 3) - 1.90) < 0.01

if 0:
    from pylab import *

    f = f0
    #f= f1

    s = linspace(0, 3, 16)
    t = '-'
    for xc in [pbe_1, pw91_1]:
        for rs in [0.0001, 2, 10]:
            plot(s, [f(xc, rs, x) for x in s], t)
        t = 'o'
    show()
