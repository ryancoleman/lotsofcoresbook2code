import os
from math import exp, pi, sqrt
import numpy as np

from gpaw.gauss import Gauss, Lorentz
from gpaw.test import equal
from gpaw.utilities.folder import Folder

# Gauss and Lorentz functions

width = 0.5
x = 1.5

equal(Gauss(width).get(x), 
      exp(- x**2 / 2 / width**2) / sqrt(2 * pi) / width, 
      1.e-15)
equal(Lorentz(width).get(x), 
      width / (x**2 + width**2) / pi, 
      1.e-15)

# folder function

for name in ['Gauss', 'Lorentz']:
    folder = Folder(width, name)

    x = [0, 2]
    y = [[2, 0, 1], [1, 1, 1]]

    xl, yl = folder.fold(x, y, dx=.7)

    # check first value
    if name == 'Lorentz':
        func = Lorentz(width)
    else:
        func = Gauss(width)
    yy = np.dot(np.array(y)[:, 0], func.get(xl[0] - np.array(x)))
    equal(yl[0, 0], yy, 1.e-15)

# write spectrum

from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft.spectrum import spectrum

fname = 'lr.dat.gz'
if os.path.exists(fname):
    lr = LrTDDFT(fname)
    lr.diagonalize()
    spectrum(lr, 'spectrum.dat')

