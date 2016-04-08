from gpaw.tddft import TDDFT
from matplotlib import use, patches
from matplotlib.pyplot import *
from ase.units import Bohr
import numpy as np
import sys


# Initialize TDDFT and QSFDTD
td_calc = TDDFT('gs.gpw')

def generate_xygrid(d, g, box):
    
    vslice = 2 # yx
    
    # Determine the array lengths in each dimension
    ng = d.shape
    
    X = None
    Y = None
    U = None
    V = None
    
    # Slice data
    d_slice=np.rollaxis(d, vslice)[g[vslice], :, :]
    d_proj = np.zeros(d_slice.shape)
    for ind, val in np.ndenumerate(d_slice):
        d_proj[ind] = np.where(np.append(np.rollaxis(d, vslice)[:, ind[0], ind[1]], 1.0)!=0)[0][0]

    # Grids
    x = np.linspace(0, box[1], ng[1])
    y = np.linspace(0, box[0], ng[0])
    
    # Meshgrid and corresponding data
    X, Y = np.meshgrid(x, y)
    U = np.real(d_slice[1]) #y
    V = np.real(d_slice[0]) #x

    # Spacing
    dx = x[1]-x[0]
    dy = y[1]-y[0]    
    return d_slice, d_proj, (x, y, dx, dy), (X, Y, U, V) 

use('Agg')

poisson_solver = td_calc.hamiltonian.poisson
atoms = td_calc.atoms

box = np.diagonal(poisson_solver.cl.gd.cell_cv) * Bohr  # in Ang
atom_positions = atoms.get_positions() + poisson_solver.qm.corner1 * Bohr  # in Ang
atom_elements  = atoms.get_chemical_symbols()

# create figure
figure(1, figsize=(4, 4))
rcParams['font.size'] = 14

# prepare data
plotData = poisson_solver.classical_material.beta[0]
ng = plotData.shape
colorLimits = None #[-1e-3, 1e-3]

axis = 2
ax = subplot(1, 1, 1)
g = [None, None, ng[2]/2]

dmy1, d_proj, (x, y, dx, dy), dmy2 = generate_xygrid(plotData, g, box)

imshow(d_proj, interpolation='bicubic', origin='lower',
      cmap=cm.Blues,
      extent=[x[0]-dx/2, x[-1]+dx/2, y[0]-dy/2, y[-1]+dy/2])

# Plot atoms
flt = np.array([True,True,False])
for position in atom_positions:
    scatter(position[1], position[0], s=50, c='r', marker='o')

xlim(x[0], x[-1])
ylim(y[0], y[-1])

# Mark the quantum region
i, j = 1, 0
qmrect = patches.Rectangle((poisson_solver.qm.corner1[i] * Bohr, poisson_solver.qm.corner1[j] * Bohr),
                           (poisson_solver.qm.corner2[i] - poisson_solver.qm.corner1[i]) * Bohr,
                           (poisson_solver.qm.corner2[j] - poisson_solver.qm.corner1[j]) * Bohr,
                           color='black',  # #0099FF',
                           fill=0,
                           linewidth=1.0)
ax.add_patch(qmrect)

# Classical
dmy1, dmy_proj, (x, y, dx, dy), dmy3 = generate_xygrid(plotData, g, box)
xx, yy = np.meshgrid(x, y)
scatter(xx,
        yy,
        s=0.75, c='k', marker='o')

# Quantum
dmy1, dmy_proj, (x, y, dx, dy), dmy3 = generate_xygrid(plotData, g, box = np.diagonal(poisson_solver.qm.gd.cell_cv) * Bohr)
xx, yy = np.meshgrid(x, y)
scatter(poisson_solver.qm.corner1[i] * Bohr + xx,
        poisson_solver.qm.corner1[j] * Bohr + yy,
        s=0.25, c='k', marker='o')

# Labels
xlabel('y [Ang]')
ylabel('x [Ang]')

# Plot
tight_layout()
savefig('geom.png')

import os
os.system('cp geom.png ../../_build')

