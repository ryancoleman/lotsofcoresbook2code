import numpy as np
from math import pi
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('figure', figsize=(4.0, 4.0), dpi=800)

data = np.loadtxt('df.csv', delimiter=',')

# Set plotting range
xmin = 0.1
xmax = 10.0
inds_w = (data[:, 0] >= xmin) & (data[:, 0] <= xmax)

plt.plot(data[inds_w, 0], 4 * pi * data[inds_w, 4])
plt.xlabel('$\\omega / [eV]$')
plt.ylabel('$\\mathrm{Im} \\epsilon(\\omega)$')
plt.savefig('si_abs.png', bbox_inches='tight')
