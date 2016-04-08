from __future__ import print_function
import pickle
import numpy as np
from gpaw import GPAW, setup_paths
from gpaw.vdw import VDWFunctional

setup_paths.insert(0, '.')

d = np.linspace(3.0, 5.5, 11)
for symbol in ['Ar', 'Kr']:
    vdw = VDWFunctional(verbose=1)
    e = np.empty(11)
    de = np.empty(11)
    for i, r in enumerate(d):
        calc = GPAW('%s-dimer-%.2f.gpw' % (symbol, r), txt=None)
        e[i] = calc.get_atoms().get_potential_energy()
        de[i] = calc.get_xc_difference(vdw)
        print(i, e, de)
    calc = GPAW('%s-atom.gpw' % symbol, txt=None)
    e0 = calc.get_atoms().get_potential_energy()
    de0 = calc.get_xc_difference(vdw)
    print(e, de, e0, de0)
    pickle.dump((e, de, e0, de0), open(symbol + '.new.pckl', 'w'))

e = np.empty(11)
de = np.empty(11)
vdw = VDWFunctional(verbose=1)
for i, r in enumerate(d):
    calc = GPAW('benzene-dimer-%.2f.gpw' % r, txt=None)
    e[i] = calc.get_atoms().get_potential_energy()
    de[i] = calc.get_xc_difference(vdw)
calc = GPAW('benzene.gpw', txt=None)
e0 = calc.get_atoms().get_potential_energy()
de0 = calc.get_xc_difference(vdw)
print(e, de, e0, de0)
pickle.dump((e, de, e0, de0), open('benzene.new.pckl', 'w'))
