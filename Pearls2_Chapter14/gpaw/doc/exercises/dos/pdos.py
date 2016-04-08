import matplotlib.pyplot as plt
from gpaw import GPAW

calc = GPAW('ferro.gpw', txt=None)

ef = calc.get_fermi_level()

# Plot s, p, d projected LDOS:
for c in 'spd':
    energies, ldos = calc.get_orbital_ldos(a=0, spin=0, angular=c, width=0.4)
    plt.plot(energies - ef, ldos, label=c + '-up')

    energies, ldos = calc.get_orbital_ldos(a=0, spin=1, angular=c, width=0.4)
    plt.plot(energies - ef, ldos, label=c + '-down')

plt.legend()
plt.show()
