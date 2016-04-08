import pickle
import numpy as np
from ase.dft.kpoints import ibz_points, get_bandpath
from gpaw import GPAW
points = ibz_points['fcc']
G = points['Gamma']
X = points['X']
W = points['W']
K = points['K']
L = points['L']
calc = GPAW('Si-PBE.gpw',
            txt=None,
            parallel={'domain': 1},
            fixdensity=True,
            symmetry='off',
            convergence={'bands': 4})
kpts, x, X = get_bandpath([W, L, G, X, W, K], calc.atoms.cell)
calc.set(kpts=kpts)
calc.get_potential_energy()
e_kn = np.array([calc.get_eigenvalues(k) for k in range(len(kpts))])
pickle.dump((x, X, e_kn), open('eigenvalues.pckl', 'w'))
