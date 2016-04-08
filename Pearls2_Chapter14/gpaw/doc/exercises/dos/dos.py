import sys
import matplotlib.pyplot as plt
from gpaw import GPAW

# The following five lines read a file name and an optional width
# from the command line.
filename = sys.argv[1]
if len(sys.argv) > 2:
    width = float(sys.argv[2])
else:
    width = None

calc = GPAW(filename, txt=None)
try:
    ef = calc.get_fermi_level()
except ValueError:
    ef = 0
energy, dos = calc.get_dos(spin=0, width=width)
plt.plot(energy - ef, dos)
if calc.get_number_of_spins() == 2:
    energy, dos = calc.get_dos(spin=1, width=width)
    plt.plot(energy - ef, dos)
    plt.legend(('up', 'down'), loc='upper left')
plt.xlabel(r'$\epsilon - \epsilon_F \ \rm{(eV)}$')
plt.ylabel('Density of States (1/eV)')
plt.show()
