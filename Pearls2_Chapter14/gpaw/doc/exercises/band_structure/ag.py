import numpy as np
from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from gpaw import GPAW

# Perform standard ground state calculation (with plane wave basis)
ag = bulk('Ag')
calc = GPAW(mode='pw',
            xc='GLLBSC',
            kpts=(10, 10, 10),
            txt='Ag_GLLBSC.txt')
ag.set_calculator(calc)
ag.get_potential_energy()
calc.write('Ag_GLLBSC.gpw')
ef = calc.get_fermi_level()

# Restart from ground state and fix potential:
calc = GPAW('Ag_GLLBSC.gpw',
            nbands=16,
            basis='dzp',
            fixdensity=True,
            symmetry='off',
            convergence={'bands': 12})

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
emax = e_kn[:, 12].max() + 1.0

plt.figure(figsize=(5, 8))
for n in range(12):
    plt.plot(x, e_kn[:, n])
for p in X:
    plt.plot([p, p], [emin, emax], 'k-')
plt.plot([0, X[-1]], [0, 0], 'k-')
plt.xticks(X, ['$%s$' % n for n in ['W', 'L', r'\Gamma', 'X', 'W', 'K']])
plt.axis(xmin=0, xmax=X[-1], ymin=emin, ymax=emax)
plt.xlabel('k-vector')
plt.ylabel('E - E$_F$ (eV)')
plt.title('Bandstructure of Silver')
plt.savefig('Ag.png')
