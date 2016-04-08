from __future__ import print_function
import numpy as np
from ase import Atom, Atoms
from ase.lattice import bulk
from ase.units import Hartree, Bohr
from gpaw import GPAW, FermiDirac
from gpaw.eigensolvers.rmm_diis_old import RMM_DIIS
from gpaw.mixer import Mixer
from gpaw.response.bse import BSE

GS = 1
bse = 1
df = 1
check_spectrum = 1

if GS:
    a = 6.75 * Bohr
    atoms = bulk('C', 'diamond', a=a)

    calc = GPAW(h=0.2,
                eigensolver=RMM_DIIS(),
                mixer=Mixer(0.1,3),
                kpts=(2,2,2),
                occupations=FermiDirac(0.001),
                nbands=8,
                convergence={'bands':'all'})

    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    calc.write('C_kpt8.gpw','all')

if bse:
    bse = BSE('C_kpt8.gpw',
              w=np.linspace(0,20,201),
              mode='RPA',
              nc=[0,8],
              nv=[0,8],
              coupling=True,
              q=np.array([0.0001,0,0.]),
              optical_limit=True,
              ecut=50.,
              nbands=8)
    bse.get_dielectric_function('C_bse.dat')

if df:
    from gpaw.response.df0 import DF
    df = DF('C_kpt8.gpw',
            w=np.linspace(0,20,201),
            q=np.array([0.0001,0,0.]),
            optical_limit=True,
            ecut=50.,
            hilbert_trans=False)
    df.get_absorption_spectrum(filename='C.dat')


if check_spectrum:
    d = np.loadtxt('C_bse.dat')[:,2]
    Nw1 = 96
    Nw2 = 109
    if d[Nw1] > d[Nw1-1] and d[Nw1] > d[Nw1+1] and \
       d[Nw2] > d[Nw2-1] and d[Nw2] > d[Nw2+1] :
        pass
    else:
        print(d[Nw1], d[Nw2])
        raise ValueError('Absorption peak not correct ! ')

    if np.abs(d[Nw1] - 67.1975750304) > 1e-5 or \
       np.abs(d[Nw2] - 90.9851151994) > 1e-5 :
        print(d[Nw1], d[Nw2])
        raise ValueError('Please check spectrum strength ! ')

    d2 = np.loadtxt('C.dat.x')
    print(np.abs(d - d2[:200, 4]).sum())
    if np.abs(d - d2[:200, 4]).sum() > 1e-3:
        raise ValueError('Please compare two spectrum')

  
