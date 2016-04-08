from __future__ import print_function
from gpaw import GPAW, restart, FermiDirac, PoissonSolver
from ase import Atoms
from gpaw.test import equal, gen
import os
from gpaw.mpi import world

gen('Si', xcname='GLLBSC')

e = {}
niter = {}

energy_tolerance = 0.0001
niter_tolerance = 0

e_ref = {'LDA': {'restart': -5.5728768784094758},
         'GLLBSC': {'restart': -5.4458036264351}} # svnversion 5252
niter_ref = {'LDA': {'restart': 16},
             'GLLBSC': {'restart': 16}} # svnversion 5252

for xc in ['LDA','GLLBSC']:
    a = 4.23
    bulk = Atoms('Si2', cell=(a, a, a), pbc=True,
              scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]])
    calc = GPAW(h=0.25,
                nbands=8,
                poissonsolver=PoissonSolver(nn='M'),
                occupations=FermiDirac(width=0.01),
                kpts=(3, 3, 3),
                convergence={'eigenstates': 9.2e-11,
                             'bands': 8},
                xc=xc,
                eigensolver='cg')

    bulk.set_calculator(calc)
    e[xc] = {'direct': bulk.get_potential_energy()}
    niter[xc] = {'direct': calc.get_number_of_iterations()}
    print(calc.get_ibz_k_points())
    old_eigs = calc.get_eigenvalues(kpt=3)
    calc.write('Si_gs.gpw')
    del bulk
    del calc
    bulk, calc = restart('Si_gs.gpw',
                         fixdensity=True,
                         kpts=[[0, 0, 0], [1.0 / 3, 1.0 / 3, 1.0 / 3]])
    e[xc] = {'restart': bulk.get_potential_energy()}
    niter[xc] = {'restart': calc.get_number_of_iterations()}

    if world.rank == 0:
        os.remove('Si_gs.gpw')
    diff = calc.get_eigenvalues(kpt=1)[:6]-old_eigs[:6]
    if world.rank == 0:
        print("occ. eig. diff.", diff)
        error = max(abs(diff))
        assert error < 5e-6

    for mode in e[xc].keys():
        equal(e[xc][mode], e_ref[xc][mode], energy_tolerance)
