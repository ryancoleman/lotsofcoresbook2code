#!/usr/bin/env python

from ase.visualize import view
from ase.lattice.surface import fcc111, add_adsorbate
from gpaw import GPAW
from gpaw.mixer import MixerSum
from gpaw import dscf

filename='homo'

#-------------------------------------------

c_mol = GPAW(nbands=9, h=0.2, xc='RPBE', kpts=(8,6,1),
             spinpol=True,
             convergence={'energy': 100,
                          'density': 100,
                          'eigenstates': 1.0e-9,
                          'bands': 'occupied'}, txt='CO_homo.txt')

calc = GPAW(nbands=45, h=0.2, xc='RPBE', kpts=(8,6,1),
            eigensolver='cg',
            spinpol=True,
            mixer=MixerSum(nmaxold=5, beta=0.1, weight=100),
            convergence={'energy': 100,
                         'density': 100,
                         'eigenstates': 1.0e-7,
                         'bands': -10}, txt=filename+'.txt')

#----------------------------------------

#  Import Slab with relaxed CO
#slab = ('gs.gpw').get_atoms()
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

#Homo wavefunction
wf_u = [kpt.psit_nG[4] for kpt in c_mol.wfs.kpt_u]

#Homo projector overlaps
mol = range(len(slab))[-2:]
p_uai = [dict([(mol[a], P_ni[4]) for a, P_ni in kpt.P_ani.items()])
         for kpt in c_mol.wfs.kpt_u]

#   Slab with adsorbed molecule
#-----------------------------------
slab.set_calculator(calc)
orbital = dscf.AEOrbital(calc, wf_u, p_uai, Estart=-100.0, Eend=0.0)
dscf.dscf_calculation(calc, [[-1.0, orbital, 1]], slab)
slab.get_potential_energy()

