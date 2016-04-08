from __future__ import print_function
from ase.io import read
from gpaw import GPAW
from gpaw import dscf
from gpaw.mixer import MixerSum
from numpy import reshape, dot
from gpaw.test import equal

filename='excited'

#-------------------------------------------

c_mol = GPAW(nbands=8, h=0.2, xc='RPBE', kpts=(4,6,1),
             eigensolver='cg', spinpol=True,
             convergence={'bands': -2}, txt='CO.txt')

calc = GPAW(nbands=120, h=0.2, xc='RPBE', kpts=(4,6,1),
            setups={'Pt': '10'},
            eigensolver='cg', spinpol=True,
            mixer=MixerSum(nmaxold=5, beta=0.1, weight=100),
            convergence={'eigenstates': 1.0e-4,
                         'bands': -10}, txt=filename+'.txt')

#----------------------------------------

#Import Slab with relaxed CO
slab = read('CO_Pt.traj')
pos = slab.get_positions()
slab.set_calculator(calc)
E_gs = slab.get_potential_energy()

#Molecule
molecule = slab.copy()
del molecule[:-2]
molecule.set_calculator(c_mol)

del slab[-2:]

d_0 = pos[-1][2] - pos[-2][2]
z_0 = 16*pos[-1][2]/(12+16) + 12*pos[-2][2]/(12+16)

ds = [d_0 - 0.02 + i*0.02 for i in range(3)]

E_es = []

for d in ds:
    molecule[0].position = (pos[12][0], pos[12][1], z_0 - 16*d/(12+16))
    molecule[1].position = (pos[13][0], pos[13][1], z_0 + 12*d/(12+16))
    molecule.get_potential_energy()

    #Find band corresponding to lumo
    lumo = c_mol.get_pseudo_wave_function(band=6, kpt=0, spin=1)
    lumo = reshape(lumo, -1)

    wf1_k = [c_mol.get_pseudo_wave_function(band=5, kpt=k, spin=1)
             for k in range(len(c_mol.wfs.kd.weight_k))]
    wf2_k = [c_mol.get_pseudo_wave_function(band=6, kpt=k, spin=1)
             for k in range(len(c_mol.wfs.kd.weight_k))]

    band_k = []
    for k in range(len(c_mol.wfs.kd.weight_k)):
        wf1 = reshape(wf1_k[k], -1)
        wf2 = reshape(wf2_k[k], -1)
        p1 = abs(dot(wf1, lumo))
        p2 = abs(dot(wf2, lumo))
        if p1 > p2:
            band_k.append(5)
        else:
            band_k.append(6)
        print('Kpt', k, p1, p2, 'band', band_k[-1])

    #Lumo wavefunction
    wf_u = [kpt.psit_nG[band_k[kpt.k]] for kpt in c_mol.wfs.kpt_u]

    #Lumo projector overlaps
    mol = [12,13]
    p_uai = [dict([(mol[a], P_ni[band_k[kpt.k]])
                   for a, P_ni in kpt.P_ani.items()])
             for kpt in c_mol.wfs.kpt_u]

    #   Slab with adsorbed molecule
    #-----------------------------------
    slab.extend(molecule)
    orbital = dscf.AEOrbital(calc, wf_u, p_uai)
    dscf.dscf_calculation(calc, [[1.0, orbital, 1]], slab)
    E_es.append(slab.get_potential_energy())

    del slab[12:]

F_d = (E_es[0]-E_es[2]) / 0.02 / 2
equal(E_es[1], E_gs + 3.9, 0.1)
equal(F_d, 3.8, 0.2)
