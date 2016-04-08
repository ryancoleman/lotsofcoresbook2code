from __future__ import print_function
from ase.lattice import bulk
from sys import argv
from ase.dft.kpoints import ibz_points, get_bandpath
from gpaw.eigensolvers.davidson import Davidson
from gpaw import *
from ase import *
from gpaw.test import gen
from gpaw import setup_paths

import os

"""This calculation has the following structure.
1) Calculate the ground state of Diamond.
2) Calculate the band structure of diamond in order to obtain accurate KS band gap for Diamond.
3) Calculate ground state again, and calculate the potential discontinuity using accurate band gap.
4) Calculate band structure again, and apply the discontinuity to CBM.

Compare to reference.
"""

xc = 'GLLBSC'
gen('C',xcname=xc)
setup_paths.insert(0, '.')

# Calculate ground state
atoms = bulk('C', 'diamond', a=3.567)
calc = GPAW(h=0.15, kpts=(4,4,4), xc=xc, nbands = 6,
            eigensolver=Davidson(niter=2))
atoms.set_calculator(calc)
atoms.get_potential_energy()
calc.write('Cgs.gpw')

# Calculate accurate KS-band gap from band structure
points = ibz_points['fcc']

# CMB is in G-X
G = points['Gamma']
X = points['X']

#W = points['W']
#K = points['K']
#L = points['L']
#[W, L, G, X, W, K]

kpts, x, X = get_bandpath([G, X], atoms.cell, npoints=12)
calc = GPAW('Cgs.gpw', kpts=kpts, fixdensity=True, symmetry='off',
            convergence=dict(bands=6), eigensolver=Davidson(niter=2))
calc.get_atoms().get_potential_energy()
# Get the accurate KS-band gap
homolumo = calc.occupations.get_homo_lumo(calc.wfs)
homo, lumo = homolumo
print("band gap ",(lumo-homo)*27.2)
    
# Redo the ground state calculation
calc = GPAW(h=0.15, kpts=(4,4,4), xc=xc, nbands = 6,
            eigensolver=Davidson(niter=2))
atoms.set_calculator(calc)
atoms.get_potential_energy()
# And calculate the discontinuity potential with accurate band gap
response = calc.hamiltonian.xc.xcs['RESPONSE']
response.calculate_delta_xc(homolumo=homolumo)
calc.write('CGLLBSC.gpw')

# Redo the band structure calculation
atoms, calc = restart('CGLLBSC.gpw', kpts=kpts, fixdensity=True, symmetry='off',
                      convergence=dict(bands=6), eigensolver=Davidson(niter=2))
atoms.get_potential_energy()
response = calc.hamiltonian.xc.xcs['RESPONSE']
KS, dxc = response.calculate_delta_xc_perturbation()

assert abs(KS+dxc-5.41)<0.10
#M. Kuisma et. al, Phys. Rev. B 82, 115106, QP gap for C, 5.41eV, expt. 5.48eV
