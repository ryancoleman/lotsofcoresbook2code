"""Test read/write of restart files between fd and lcao mode"""
from __future__ import print_function
import os

from ase import Atom, Atoms
from gpaw import GPAW, Mixer, restart, FermiDirac
from gpaw.test import equal

energy_tolerance = 0.0001
niter_tolerance = 0

if not os.path.isfile('Na4_fd.gpw'):
    # Do grid kpts calculation
    a = 3.31
    atoms = Atoms([Atom('Na',(i*a,0,0)) for i in range(4)], pbc=(1,0,0))
    atoms.center(vacuum=a/2, axis=0)
    atoms.center(vacuum=3.5, axis=1)
    atoms.center(vacuum=3.5, axis=2)

    calc = GPAW(nbands=-3,
                h=0.3,
                setups={'Na': '1'},
                xc='PBE',
                occupations=FermiDirac(width=0.1),
                kpts=(3, 1, 1),
                #basis='dzp',
                txt='Na4_fd.txt')
    atoms.set_calculator(calc)
    etot_fd = atoms.get_potential_energy()
    niter_fd = calc.get_number_of_iterations()
    print('Etot:', etot_fd, 'eV in fd-mode')
    calc.write('Na4_fd.gpw')
    del atoms,calc

    equal(etot_fd, -1.99055, energy_tolerance)

if os.path.isfile('Na4_fd.gpw'):
    # LCAO calculation based on grid kpts calculation
    atoms, calc = restart('Na4_fd.gpw',
                          #basis='dzp',
                          mode='lcao',
                          txt='Na4_lcao.txt')
    etot_lcao = atoms.get_potential_energy()
    niter_lcao = calc.get_number_of_iterations()
    print('Etot:', etot_lcao, 'eV in lcao-mode')
    calc.write('Na4_lcao.gpw')
    del atoms, calc

    equal(etot_lcao, -1.9616, energy_tolerance)
