import numpy as np
import matplotlib.pyplot as plt

from gpaw.response.chi0 import frequency_grid


omegamax = 50.0
domega0 = 0.2

plt.figure(figsize=(5, 5))
for omega2 in [2.5, 5, 10, 15, 20, np.inf]:
    x = frequency_grid(domega0, omega2, omegamax)
    y = range(len(x))
    if omega2 == np.inf:
        label = '$\\omega_2 = \\infty$'
    else:
        label = '$\\omega_2 = %.1f\\, \\mathrm{eV}$' % omega2
    plt.plot(x, y, '.', label=label)
plt.ylabel('Freq. no')
plt.xlabel('$\\omega\\, [\mathrm{eV}]$')
plt.axis(xmax=30, ymax=200)
plt.title('$\\Delta\\omega_0 = 0.2\\, \mathrm{eV}$')
plt.legend(loc=2)
plt.savefig('nl_freq_grid.png', bbox_inches='tight')
