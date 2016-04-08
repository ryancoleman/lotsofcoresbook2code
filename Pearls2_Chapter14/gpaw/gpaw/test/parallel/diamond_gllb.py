from __future__ import print_function
from ase.lattice import bulk
from sys import argv
from ase.dft.kpoints import ibz_points, get_bandpath
from gpaw import *
from ase import *
from gpaw.test import gen
from gpaw import setup_paths

import os

"""
Calculate diamond with various parallelizations with GLLBSC
"""

xc = 'GLLBSC'
gen('C',xcname=xc)
setup_paths.insert(0, '.')

KSb = []
dxcb = []

eigensolver='rmm-diis'

for band in [1,2,4]:
    # Calculate ground state
    atoms = bulk('C', 'diamond', a=3.567)
    calc = GPAW(h=0.15, kpts=(4,4,4), xc=xc, nbands = 8,
                eigensolver=eigensolver,
                parallel={'band':band})
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    calc.write('Cgs.gpw')

    # Calculate accurate KS-band gap from band structure
    points = ibz_points['fcc']

    # CMB is in G-X
    G = points['Gamma']
    X = points['X']

    kpts, x, X = get_bandpath([G, X], atoms.cell, npoints=12)
    calc = GPAW('Cgs.gpw', kpts=kpts, fixdensity=True, symmetry='off',
                nbands=8, convergence=dict(bands=8), 
                eigensolver=eigensolver,
                parallel={'band':band})
    calc.get_atoms().get_potential_energy()
    # Get the accurate KS-band gap
    homolumo = calc.occupations.get_homo_lumo(calc.wfs)
    homo, lumo = homolumo
    print("band gap ",(lumo-homo)*27.2)
    
    # Redo the ground state calculation
    calc = GPAW(h=0.15, kpts=(4,4,4), xc=xc, nbands = 8, 
                eigensolver=eigensolver,
                parallel={'band':band})
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    # And calculate the discontinuity potential with accurate band gap
    response = calc.hamiltonian.xc.xcs['RESPONSE']
    response.calculate_delta_xc(homolumo=homolumo)
    calc.write('CGLLBSC.gpw')

    # Redo the band structure calculation
    atoms, calc = restart('CGLLBSC.gpw', kpts=kpts, fixdensity=True,
                          symmetry='off', convergence=dict(bands=8))
    atoms.get_potential_energy()
    response = calc.hamiltonian.xc.xcs['RESPONSE']
    KS, dxc = response.calculate_delta_xc_perturbation()
    KSb.append(KS)
    dxcb.append(dxc)    
    assert abs(KS+dxc-5.41)<0.10
    #M. Kuisma et. al, Phys. Rev. B 82, 115106, QP gap for C, 5.41eV, expt. 5.48eV

assert abs(KSb[0]-KSb[1])<1e-6
assert abs(KSb[0]-KSb[2])<1e-6
assert abs(dxcb[0]-dxcb[1])<1e-6
assert abs(dxcb[0]-dxcb[2])<1e-6
