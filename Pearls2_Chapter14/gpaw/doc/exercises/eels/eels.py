from gpaw import GPAW
from gpaw.response.df import DielectricFunction

calc = GPAW('Ag_GLLBSC.gpw')
calc.diagonalize_full_hamiltonian(nbands=30)
calc.write('Ag_GLLBSC_full.gpw', 'all')

# Set up dielectric function:
df = DielectricFunction(calc='Ag_GLLBSC_full.gpw',  # Ground state input
                        domega0=0.05)  # energy grid spacing at omega=0

# Momentum transfer, must be the difference between two kpoints!
q_c = [1.0 / 10, 0, 0]
df.get_eels_spectrum(q_c=q_c)  # a file called 'eels.csv' is generated

# Plot spectrum
import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('eels.csv', delimiter=',')
omega = data[:, 0]
eels = data[:, 2]
plt.plot(omega, eels)
plt.xlabel('Energy (eV)')
plt.ylabel('Loss spectrum')
plt.xlim(0, 20)
plt.show()
