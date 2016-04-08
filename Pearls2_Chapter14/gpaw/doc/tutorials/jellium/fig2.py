import numpy as np
# mathtext fails to create title with matplotlib 0.99 on el6
import matplotlib
matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as plt
from ase.units import Bohr
from gpaw import GPAW
rs = 5.0 * Bohr
calc = GPAW('surface.gpw', txt=None)
density = calc.get_pseudo_density()[0, 0]
h = 0.2
a = 8 * h
v = 3 * a
L = 10 * a
z = np.linspace(0, v + L + v, len(density), endpoint=False)
# Position of surface is between two grid points:
z0 = (v + L - h / 2)
n = 1 / (4 * np.pi / 3 * rs**3)  # electron density
kF = (3 * np.pi**2 * n)**(1.0 / 3)
lambdaF = 2 * np.pi / kF  # Fermi wavelength
plt.figure(figsize=(6, 6 / 2**0.5))
plt.plot([-L / lambdaF, -L / lambdaF, 0, 0], [0, 1, 1, 0], 'k')
plt.plot((z - z0) / lambdaF, density / n)
#plt.xlim(xmin=-1.2, xmax=1)
plt.ylim(ymin=0)
plt.title(r'$r_s=%.1f\ \mathrm{Bohr},\ \lambda_F=%.1f\ \mathrm{Bohr}$' %
          (rs / Bohr, lambdaF / Bohr))
plt.xlabel('DISTANCE (FERMI WAVELENGTHS)')
plt.ylabel('ELECTRON DENSITY')
plt.savefig('fig2.png')

