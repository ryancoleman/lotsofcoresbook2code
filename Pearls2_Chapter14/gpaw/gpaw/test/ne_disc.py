import os
from ase import *
from gpaw import GPAW, restart
from gpaw.atom.generator import Generator
from gpaw.atom.configurations import parameters
from gpaw import setup_paths
from gpaw.mpi import world
from gpaw.test import equal

# This test calculates the derivative discontinuity of Ne-atom
# first on 3D without restart. Then does restart and recalculates.

atom = 'Ne'
setup_paths.insert(0, '.')

for xcname in ['GLLB','GLLBSC']:
    if world.rank == 0:
        g = Generator(atom, xcname =xcname, scalarrel=False,nofiles=True)
        g.run(**parameters[atom])
        eps = g.e_j[-1]
    world.barrier()

    a = 10
    Ne = Atoms([Atom(atom, (0, 0, 0))],
               cell=(a, a, a), pbc=False)
    Ne.center()
    calc = GPAW(nbands=10, h=0.21, xc=xcname)
    Ne.set_calculator(calc)
    e = Ne.get_potential_energy()
    response = calc.hamiltonian.xc.xcs['RESPONSE']
    response.calculate_delta_xc()
    KS, dxc = response.calculate_delta_xc_perturbation()
    if xcname=='GLLB':
        equal(KS+dxc, 24.71, 1.5e-1)
    else:
        equal(KS+dxc, 27.70, 6.0e-2)
    eps3d = calc.wfs.kpt_u[0].eps_n[3]
    calc.write('Ne_temp.gpw')

    atoms, calc = restart('Ne_temp.gpw')
    KS2, dxc2 = response.calculate_delta_xc_perturbation()
    equal(KS, KS2, 1e-5)
    equal(dxc2, dxc, 1e-5)
    
    # Hardness of Ne 24.71eV by GLLB+Dxc, experimental I-A = I = 21.56eV

    if world.rank == 0:
        equal(eps, eps3d, 1e-3)
    if xcname=='GLLB':
        equal(24.71, KS2+dxc2, 1.2e-1)
