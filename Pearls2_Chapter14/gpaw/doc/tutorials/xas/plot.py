from gpaw import GPAW, setup_paths
from gpaw.xas import XAS
import pylab as plt

setup_paths.insert(0, '.')

dks_energy = 532.774  # from dks calcualtion

calc = GPAW('h2o_xas.gpw')

xas = XAS(calc, mode='xas')
x, y = xas.get_spectra(fwhm=0.5, linbroad=[4.5, -1.0, 5.0])
x_s, y_s = xas.get_spectra(stick=True)

shift = dks_energy - x_s[0]  # shift the first transition 

y_tot = y[0] + y[1] + y[2]
y_tot_s = y_s[0] + y_s[1] + y_s[2]

plt.plot(x + shift, y_tot)
plt.bar(x_s + shift, y_tot_s, width=0.001)
plt.show()
