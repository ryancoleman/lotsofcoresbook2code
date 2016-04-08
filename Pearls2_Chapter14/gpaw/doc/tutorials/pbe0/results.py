from __future__ import print_function
import pickle
import numpy as np
import matplotlib.pyplot as plt
from ase.utils.eos import EquationOfState
import ase.units as units

results = []
K = range(3, 13)
for k in K:
    A, data = pickle.load(open('eos-%d.pckl' % k))
    for energies, xcname in zip(data, ['PBE', 'LDA', 'PBE0']):
        eos = EquationOfState(A**3 / 4, energies)
        v0, e0, B = eos.fit()
        a = (v0 * 4)**(1 / 3.0)
        B *= 1.0e24 / units.kJ
        print(('%-4s %2d %.3f %.3f %.2f' % (xcname, k, a, e0, B)))
        results.append((a, B))

results = np.array(results).reshape((-1, 3, 2))

LDA = dict(
    a=5.4037,
    B=95.1,
    eGX=0.52)
PBE = dict(
    a=5.469,
    B=87.8,
    eGG=2.56,
    eGX=0.71,
    eGL=1.54,
    eI=0.47, # indirect
    ea=4.556)
PBE0 = dict(
    a=5.433,
    B=99.0,
    eGG=3.96,
    eGX=1.93,
    eGL=2.87,
    eI=1.74,
    ea=4.555)

plt.figure(figsize=(6, 4))
plt.plot(K, results[:, 0, 0], label='PBE')
plt.plot(K, results[:, 1, 0], label='LDA')
plt.plot(K, results[:, 2, 0], label='PBE0')
plt.plot([12], [PBE['a']], 'b*', label='_nolegend_')
plt.plot([12], [LDA['a']], 'g*', label='_nolegend_')
plt.plot([12], [PBE0['a']], 'r*', label='_nolegend_')
plt.legend(loc='upper right')
plt.xlabel('number of k-points')
plt.ylabel('lattice constant [Ang]')
plt.savefig('a.png')

plt.figure(figsize=(6, 4))
plt.plot(K, results[:, 0, 1], label='PBE')
plt.plot(K, results[:, 1, 1], label='LDA')
plt.plot(K, results[:, 2, 1], label='PBE0')
plt.plot([12], [PBE['B']], 'b*', label='_nolegend_')
plt.plot([12], [LDA['B']], 'g*', label='_nolegend_')
plt.plot([12], [PBE0['B']], 'r*', label='_nolegend_')
plt.legend(loc='upper right')
plt.xlabel('number of k-points')
plt.ylabel('bulk modulus [GPa]')
plt.savefig('B.png')

