import sys
import pickle
import numpy as np
from gpaw.xc.hybridk import HybridXC
from gpaw import GPAW
from si_pbe import groundstate
a0 = 5.43
A = np.linspace(a0 - 0.06, a0 + 0.06, 5)
k = int(sys.argv[1])
eos = np.zeros((3, 5))
for i, a in enumerate(A):
    si = groundstate(a, k)
    epbe = si.get_potential_energy()
    elda = epbe + si.calc.get_xc_difference('LDA')
    pbe0 = HybridXC('PBE0', alpha=5.0)
    epbe0 = epbe + si.calc.get_xc_difference(pbe0)
    eos[:, i] = epbe, elda, epbe0
pickle.dump((A, eos), open('eos-%d.pckl' % k, 'w'))
