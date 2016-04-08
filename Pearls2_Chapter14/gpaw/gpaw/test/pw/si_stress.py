from __future__ import print_function
import numpy as np
from ase.lattice import bulk
from gpaw import GPAW, PW


si = bulk('Si')
k = 3
si.calc = GPAW(mode=PW(250),
               xc='PBE',
               kpts=(k, k, k),
               convergence={'energy': 1e-8},
               txt='si.txt')

si.set_cell(np.dot(si.cell,
                   [[1.02, 0, 0.03],
                    [0, 0.99, -0.02],
                    [0.2, -0.01, 1.03]]),
            scale_atoms=True)

sigma_vv = si.get_stress(voigt=False)
print(sigma_vv)

deps = 1e-5
cell = si.cell.copy()
for v in range(3):
    x = np.eye(3)
    x[v, v] += deps
    si.set_cell(np.dot(cell, x), scale_atoms=True)
    ep = si.calc.get_potential_energy(si, force_consistent=True)
    x[v, v] -= 2 * deps
    si.set_cell(np.dot(cell, x), scale_atoms=True)
    em = si.calc.get_potential_energy(si, force_consistent=True)
    s = (ep - em) / 2 / deps / si.get_volume()
    print((v, s, abs(s - sigma_vv[v, v])))
    assert abs(s - sigma_vv[v, v]) < 1e-4
for v1 in range(3):
    v2 = (v1 + 1) % 3
    x = np.eye(3)
    x[v1, v2] = deps
    x[v2, v1] = deps
    si.set_cell(np.dot(cell, x), scale_atoms=True)
    ep = si.calc.get_potential_energy(si, force_consistent=True)
    x[v1, v2] = -deps
    x[v2, v1] = -deps
    si.set_cell(np.dot(cell, x), scale_atoms=True)
    em = si.calc.get_potential_energy(si, force_consistent=True)
    s = (ep - em) / deps / 4 / si.get_volume()
    print((v1, v2, s, abs(s - sigma_vv[v1, v2])))
    assert abs(s - sigma_vv[v1, v2]) < 2e-4
