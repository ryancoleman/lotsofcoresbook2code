import numpy as np
from pylab import *

A = loadtxt('frequency_equidistant.dat').transpose()
plot(A[0], A[1], label='Equidistant')
B = loadtxt('frequency_gauss16.dat').transpose()
plot(B[0], B[1], 'o', label='Gauss-Legendre 16')

xlabel('Frequency [eV]', fontsize=18)
ylabel('Integrand', fontsize=18)
axis([0, 50, None, None])
legend(loc='lower right')
#show()
savefig('E_w.png')
