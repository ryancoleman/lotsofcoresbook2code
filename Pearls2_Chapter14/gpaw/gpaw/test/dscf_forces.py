from __future__ import print_function
import sys

import numpy as np

from ase import Atom, Atoms
from gpaw import GPAW, MixerSum, dscf
from gpaw.mpi import rank

L = 10.0

eigensolver = 'cg'
xc = 'PBE'

E0 = []
E1 = []
d0 = 1.4
delta = 0.05
for i in range(3):
    d = d0 + (i - 1) * delta
    atoms = Atoms([Atom('N',[0.5*L, 0.5*L, 0.5*L - 0.5*d]),
                   Atom('N',[0.5*L, 0.5*L, 0.5*L + 0.5*d])])
    atoms.set_cell([L,L,L])

    calc_gs = GPAW(h = 0.2, 
                   nbands = -5,
                   xc = xc, 
                   width = 0.05,
                   eigensolver = eigensolver,
                   spinpol = True,
                   txt = None,
                   mixer = MixerSum(beta=0.1, nmaxold=5, weight=50.0),
                   convergence = {'energy': 100,
                                  'density': 100,
                                  'eigenstates': 1.0e-9,
                                  'bands': -1})
    atoms.set_calculator(calc_gs)

    E0.append(atoms.get_potential_energy())
    if i==1:
        F0 = atoms.get_forces()

    calc_es = GPAW(h = 0.2, 
                   nbands = -5,
                   xc = xc, 
                   width = 0.05,
                   eigensolver = eigensolver,
                   spinpol = True,
                   txt = None,
                   mixer = MixerSum(beta=0.1, nmaxold=5, weight=50.0),
                   convergence = {'energy': 100,
                                  'density': 100,
                                  'eigenstates': 1.0e-9,
                                  'bands': -1})

    n = 5 # LUMO
    molecule = [0, 1]
    wf_u = [kpt.psit_nG[n] for kpt in calc_gs.wfs.kpt_u]
    p_uai = [dict([(molecule[a], P_ni[n]) for a, P_ni in kpt.P_ani.items()])
             for kpt in calc_gs.wfs.kpt_u]

    atoms.set_calculator(calc_es)
    lumo = dscf.AEOrbital(calc_es, wf_u, p_uai)
    dscf.dscf_calculation(calc_es, [[1.0, lumo, 1]], atoms)

    E1.append(atoms.get_potential_energy())
    if i==1:
        F1 = atoms.get_forces()

f0 = np.sqrt(((F0[1,:] - F0[0,:])**2).sum()) * 0.5
f0b = (E0[1] - E0[0]) / delta     # backward gradient
f0f = (E0[2] - E0[1]) / delta     # forward gradient
if rank == 0:
    print('Ground state')
    print(E0)
    print(f0b, '<', f0, '<', f0f)
assert f0 > f0b
assert f0 < f0f 

f1 = np.sqrt(((F1[1,:] - F1[0,:])**2).sum()) * 0.5
f1b = (E1[1] - E1[0]) / delta     # backward gradient
f1f = (E1[2] - E1[1]) / delta     # forward gradient 
if rank == 0:
    print('Excited state')
    print(E1)
    print(f1b, '<', f1, '<', f1f)
assert f1 > f1b
assert f1 < f1f
