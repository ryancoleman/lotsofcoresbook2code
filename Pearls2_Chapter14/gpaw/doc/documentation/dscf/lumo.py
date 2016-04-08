#!/usr/bin/env python

from numpy import reshape, dot

from ase.visualize import view
from ase.lattice.surface import fcc111, add_adsorbate
from gpaw import GPAW
from gpaw.mixer import MixerSum
from gpaw import dscf

filename='lumo'

#-------------------------------------------

c_mol = GPAW(nbands=9, h=0.2, xc='RPBE', kpts=(8,6,1),
             spinpol=True,
             convergence={'energy': 100,
                          'density': 100,
                          'eigenstates': 1.0e-9,
                          'bands': -2}, txt='CO_lumo.txt')

calc = GPAW(nbands=60, h=0.2, xc='RPBE', kpts=(8,6,1),
            eigensolver='cg',
            spinpol=True,
            mixer=MixerSum(nmaxold=5, beta=0.1, weight=100),
            convergence={'energy': 100,
                         'density': 100,
                         'eigenstates': 1.0e-7,
                         'bands': -10}, txt=filename+'.txt')

#----------------------------------------

#  Import Slab with relaxed CO
#slab = Calculator('gs.gpw').get_atoms()
slab = fcc111('Pt', size=(1, 2, 3), orthogonal=True)
add_adsorbate(slab, 'C', 2.0, 'ontop')
add_adsorbate(slab, 'O', 3.15, 'ontop')
slab.center(axis=2, vacuum=4.0)

view(slab)

molecule = slab.copy()

del molecule [:-2]

#   Molecule
#----------------
molecule.set_calculator(c_mol)
molecule.get_potential_energy()

#Find band corresponding to lumo
lumo = c_mol.get_pseudo_wave_function(band=5, kpt=0, spin=1)
lumo = reshape(lumo, -1)

wf1_k = [c_mol.get_pseudo_wave_function(band=5, kpt=k, spin=1)
         for k in range(len(c_mol.wfs.weight_k))]
wf2_k = [c_mol.get_pseudo_wave_function(band=6, kpt=k, spin=1)
         for k in range(len(c_mol.wfs.weight_k))]

band_k = []
for k in range(len(c_mol.wfs.weight_k)):
    
    wf1 = reshape(wf1_k[k], -1)
    wf2 = reshape(wf2_k[k], -1)
    p1 = abs(dot(wf1, lumo))
    p2 = abs(dot(wf2, lumo))
    if p1 > p2:
        band_k.append(5)
    else:
        band_k.append(6)

#Lumo wavefunction
wf_u = [kpt.psit_nG[band_k[kpt.k]] for kpt in c_mol.wfs.kpt_u]

#Lumo projector overlaps
mol = range(len(slab))[-2:]
p_uai = [dict([(mol[a], P_ni[band_k[kpt.k]]) for a, P_ni in kpt.P_ani.items()])
         for kpt in c_mol.wfs.kpt_u]

#   Slab with adsorbed molecule
#-----------------------------------
slab.set_calculator(calc)
orbital = dscf.AEOrbital(calc, wf_u, p_uai)
dscf.dscf_calculation(calc, [[1.0, orbital, 1]], slab)
slab.get_potential_energy()
