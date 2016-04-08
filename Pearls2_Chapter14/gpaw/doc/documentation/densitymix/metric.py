# creates: metric.png

import numpy as np
import matplotlib
import pylab as plt
from math import pi, cos

# Special points in the BZ of a simple cubic cell
G = pi * np.array([0., 0., 0.])
R = pi * np.array([1., 1., 1.])
X = pi * np.array([1., 0., 0.])
M = pi * np.array([1., 1., 0.])

# The path for the band plot
path = [X, G, R, X, M, G]
textpath = [r'$X$', r'$\Gamma$', r'$R$', r'$X$', r'$M$', r'$\Gamma$']

# Make band data
qvec = []
lines = [0]
previous = path[0]
for next in path[1:]:
    Npoints = int(round(20 * np.linalg.norm(next - previous)))
    lines.append(lines[-1] + Npoints)
    for t in np.linspace(0, 1, Npoints):
        qvec.append((1 - t) * previous + t * next)
    previous = next

vasp = [1 / max(np.linalg.norm(q), 1e-6)**2 for q in qvec]
gpaw = [( 1 + cos(qx) + cos(qy) + cos(qz) +
        cos(qx) * cos(qy) + cos(qx) * cos(qz) + cos(qy) * cos(qz) +
        cos(qx) * cos(qy) * cos(qz)) / 8. for qx, qy, qz in qvec]

# Plot band data
fig = plt.figure(1, figsize=(5, 3), dpi=90)
fig.subplots_adjust(left=.1, right=.95)
lim = [0, lines[-1], 0, 1.25]
plt.plot(vasp, 'k:', label='VASP')
plt.plot(gpaw, 'k-', label='GPAW')
for q in lines:
    plt.plot([q, q], lim[2:], 'k-')
plt.xticks(lines, textpath)
plt.yticks([0, 1], [r'$1$', r'$w+1$'])
plt.axis(lim)

# The pad keyword to legend was deprecated in MPL v. 0.98.4
if matplotlib.__version__ < '0.98.4':
    kwpad = {'pad': 0.1, 'axespad': 0.06}
else:
    kwpad = {'borderpad': 0.2, 'borderaxespad': 0.06}

plt.legend(loc='upper right', **kwpad)
plt.title('Special metric for density changes')
plt.savefig('metric.png', dpi=90)
#plt.show()
