import pickle
import numpy as np
from ase.dft.wannier import Wannier
from ase.units import Hartree
from ase.dft.kpoints import ibz_points, get_bandpath
from gpaw import GPAW
from gpaw.xc.tools import vxc
from gpaw.xc.hybridk import HybridXC

calc = GPAW('Si-PBE.gpw', txt=None)

pbe0 = HybridXC('PBE0', alpha=5.0)
de_skn = vxc(calc, pbe0) - vxc(calc, 'PBE')
de_kn = de_skn[0, calc.wfs.kd.bz2ibz_k]

calc.wfs.ibz2bz(calc.atoms)

w = Wannier(4, calc, 'rotations.pckl')

points = ibz_points['fcc']
G = points['Gamma']
X = points['X']
W = points['W']
K = points['K']
L = points['L']
kpts, x, X = get_bandpath([W, L, G, X, W, K], calc.atoms.cell)

epbe_kn = np.array([np.linalg.eigh(w.get_hamiltonian_kpoint(kpt))[0]
                    for kpt in kpts])

for kpt in calc.wfs.kpt_u:
    kpt.eps_n += de_kn[kpt.k] / Hartree

epbe0_kn = np.array([np.linalg.eigh(w.get_hamiltonian_kpoint(kpt))[0]
                     for kpt in kpts])

pickle.dump((x, X, epbe_kn, epbe0_kn), open('eigenvalues0.pckl', 'w'))
