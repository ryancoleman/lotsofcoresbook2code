from ase.utils.extrapolate import extrapolate
import numpy as np
from pylab import *

a = np.loadtxt('rpa_N2.dat')
ext, A, B, sigma = extrapolate(a[:,0], a[:,1], reg=3, plot=False)
plot(a[:, 0]**(-1.5), a[:, 1], 'o', label='Calculated points')
es = np.array([e for e in a[:, 0]] + [10000])
plot(es**(-1.5), A + B * es**(-1.5), '--', label='Linear regression')

t = [int(a[i, 0]) for i in range(len(a))]
xticks(a[:, 0]**(-1.5), t, fontsize=12)
axis([0., 150**(-1.5), None, -4.])
xlabel('Cutoff energy [eV]', fontsize=18)
ylabel('RPA correlation energy [eV]', fontsize=18)
legend(loc='lower right')
#show()
savefig('extrapolate.png')
