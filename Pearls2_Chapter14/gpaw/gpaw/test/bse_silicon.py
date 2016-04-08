from __future__ import print_function
import os
import numpy as np
from ase import Atom, Atoms
from ase.lattice import bulk
from ase.units import Hartree, Bohr
from gpaw import GPAW, FermiDirac
from gpaw.response.bse import BSE
from ase.dft.kpoints import monkhorst_pack
from gpaw.mpi import rank

GS = 1
bse = 1
check = 1

if GS:
    kpts = (4,4,4)
 
    a = 5.431 # From PRB 73,045112 (2006)
    atoms = bulk('Si', 'diamond', a=a)
    calc = GPAW(h=0.2,
                kpts=kpts,
                occupations=FermiDirac(0.001),
                nbands=12,
                convergence={'bands':-4})
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    calc.write('Si.gpw','all')

if bse:
    eshift = 0.8
    bse = BSE('Si.gpw',
              w=np.linspace(0,10,201),
              q=np.array([0.0001, 0, 0.0]),
              optical_limit=True,
              ecut=50.,
              nc=np.array([4,6]),
              nv=np.array([2,4]),
              eshift=eshift,
              nbands=8)
    bse.get_dielectric_function('Si_bse.dat')

    if rank == 0 and os.path.isfile('phi_qaGp'):
        os.remove('phi_qaGp')

if check:
    d = np.loadtxt('Si_bse.dat')

    Nw1 = 64
    Nw2 = 77
    if d[Nw1, 2] > d[Nw1-1, 2] and d[Nw1, 2] > d[Nw1+1, 2] \
            and d[Nw2, 2] > d[Nw2-1, 2] and d[Nw2, 2] > d[Nw2+1, 2]:
        pass
    else:
        raise ValueError('Absorption peak not correct ! ')

    if np.abs(d[Nw1, 2] - 53.3382894891) > 1. \
        or np.abs(d[Nw2, 2] - 62.7667801949 ) > 2.:
        print(d[Nw1, 2], d[Nw2, 2])
        raise ValueError('Please check spectrum strength ! ')

