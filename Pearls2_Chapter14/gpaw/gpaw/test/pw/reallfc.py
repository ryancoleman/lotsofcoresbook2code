from __future__ import print_function
import numpy as np

from gpaw.test import equal
from gpaw.grid_descriptor import GridDescriptor
from gpaw.spline import Spline
import gpaw.mpi as mpi
from gpaw.lfc import LocalizedFunctionsCollection as LFC
from gpaw.wavefunctions.pw import PWDescriptor, PWLFC
from gpaw.kpt_descriptor import KPointDescriptor


x = 2.0
rc = 3.5
r = np.linspace(0, rc, 100)

n = 40
a = 8.0
gd = GridDescriptor((n, n, n), (a, a, a), comm=mpi.serial_comm)

spos_ac = np.array([(0.15, 0.45, 0.95)])

pd = PWDescriptor(45, gd, complex)
pdr = PWDescriptor(45, gd)

from gpaw.fftw import FFTPlan
print(FFTPlan)

for l in range(4):
    print(l)
    s = Spline(l, rc, 2 * x**1.5 / np.pi * np.exp(-x * r**2))

    lfc = PWLFC([[s]], pd)
    lfcr = PWLFC([[s]], pdr)
    
    c_axi = {0: np.zeros((1, 2 * l + 1), complex)}
    c_axi[0][0, 0] = 1.9
    cr_axi = {0: np.zeros((1, 2 * l + 1))}
    cr_axi[0][0, 0] = 1.9

    b = pd.zeros(1, dtype=complex)
    br = pdr.zeros(1)

    lfc.set_positions(spos_ac)
    lfc.add(b, c_axi)
    lfcr.set_positions(spos_ac)
    lfcr.add(br, cr_axi)

    a = pd.ifft(b)
    ar = pdr.ifft(br)
    equal(abs(a-ar).max(), 0, 1e-14)

    if l == 0:
        a = a[:, ::-1].copy()
        b0 = pd.fft(a)
        br0 = pdr.fft(a.real)

    lfc.integrate(b0, c_axi)
    lfcr.integrate(br0, cr_axi)
    assert abs(c_axi[0][0]-cr_axi[0][0]).max() < 1e-14

    c_axiv = {0: np.zeros((1, 2 * l + 1, 3), complex)}
    cr_axiv = {0: np.zeros((1, 2 * l + 1, 3))}
    lfc.derivative(b0, c_axiv)
    lfcr.derivative(br0, cr_axiv)
    assert abs(c_axiv[0][0]-cr_axiv[0][0]).max() < 1e-14
