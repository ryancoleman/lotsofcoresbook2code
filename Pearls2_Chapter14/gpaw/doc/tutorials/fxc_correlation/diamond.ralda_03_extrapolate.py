from ase.utils.extrapolate import extrapolate
import numpy as np

a = np.loadtxt('diamond.ralda.rpa.dat')
b = np.loadtxt('CO.ralda_rpa_C.dat')
ext, A, B, sigma = extrapolate(a[:, 0], a[:, 1] / 2 - b[:, 1],
                               reg=3, plot=False)

a = np.loadtxt('diamond.ralda.rapbe.dat')
b = np.loadtxt('CO.ralda_rapbe_C.dat')
ext, A, B, sigma = extrapolate(a[:, 0], a[:, 1] / 2 - b[:, 1],
                               reg=3, plot=False)
