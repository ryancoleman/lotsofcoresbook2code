"""Band structure tutorial

Calculate the band structure of Si along high symmetry directions
Brillouin zone
"""

import numpy as np
from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from gpaw import GPAW, PW, FermiDirac

# Perform standard ground state calculation (with plane wave basis)
si = bulk('Si', 'diamond', 5.43)
calc = GPAW(mode=PW(200),
            xc='PBE',
            kpts=(8, 8, 8),
            occupations=FermiDirac(0.01),
            txt='Si_gs.txt')
si.set_calculator(calc)
si.get_potential_energy()
calc.write('Si_gs.gpw')
ef = calc.get_fermi_level()

# Restart from ground state and fix potential:
calc = GPAW('Si_gs.gpw',
            nbands=16,
            fixdensity=True,
            symmetry='off',
            convergence={'bands': 8})

# Use ase.dft module for obtaining k-points along high symmetry directions
points = ibz_points['fcc']
G = points['Gamma']
X = points['X']
W = points['W']
K = points['K']
L = points['L']
kpts, x, X = get_bandpath([W, L, G, X, W, K], calc.atoms.cell, npoints=60)
calc.set(kpts=kpts)
calc.get_potential_energy()
e_kn = np.array([calc.get_eigenvalues(k) for k in range(len(kpts))])

# Plot the band structure
import matplotlib.pyplot as plt

e_kn -= ef
emin = e_kn.min() - 1.0
emax = e_kn[:, 8].max() + 1.0

plt.figure(figsize=(5, 6))
for n in range(8):
    plt.plot(x, e_kn[:, n])
for p in X:
    plt.plot([p, p], [emin, emax], 'k-')
plt.plot([0, X[-1]], [0, 0], 'k-')
plt.xticks(X, ['$%s$' % n for n in ['W', 'L', r'\Gamma', 'X', 'W', 'K']])
plt.axis(xmin=0, xmax=X[-1], ymin=emin, ymax=emax)
plt.xlabel('k-vector')
plt.ylabel('E - E$_F$ (eV)')
plt.title('Bandstructure of Silicon')
plt.savefig('bandstructure.png')
