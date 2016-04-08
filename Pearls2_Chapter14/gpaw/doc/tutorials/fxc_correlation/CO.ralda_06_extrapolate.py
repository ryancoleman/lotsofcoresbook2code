from ase.utils.extrapolate import extrapolate
import numpy as np

CO_rpa = np.loadtxt('CO.ralda_rpa_CO.dat')
C_rpa = np.loadtxt('CO.ralda_rpa_C.dat')
O_rpa = np.loadtxt('CO.ralda_rpa_O.dat')

a = CO_rpa
a[:, 1] -= (C_rpa[:, 1] + O_rpa[:, 1])
ext, A, B, sigma = extrapolate(a[:, 0], a[:, 1], reg=3, plot=False)
