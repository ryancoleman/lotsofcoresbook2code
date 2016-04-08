from __future__ import print_function
import numpy as np
from ase.units import Bohr, Hartree
from ase.lattice import bulk
from gpaw import GPAW
from gpaw.eigensolvers.rmm_diis_old import RMM_DIIS
from gpaw.mixer import Mixer
from gpaw.response.df0 import DF
from gpaw.response.bse import BSE

GS = 1
df = 1
bse = 1
check_spectrum = 1

if GS:
    a = 4.043
    atoms = bulk('Al', 'fcc', a=a)
    atoms.center()
    calc = GPAW(h=0.2,
                eigensolver=RMM_DIIS(),
                mixer=Mixer(0.1,3),
                kpts=(4,2,2),
                xc='LDA',
                nbands=4,
                convergence={'bands':'all'})
    
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    calc.write('Al.gpw','all')

if bse:
    bse = BSE('Al.gpw',
              w=np.linspace(0,24,241),
              nv=[0,4],
              nc=[0,4],
              coupling=True,
              mode='RPA',
              q=np.array([0.25, 0, 0]),
              ecut=50.,
              eta=0.2)
    bse.get_dielectric_function('Al_bse.dat')
    
if df:
    # Excited state calculation
    q = np.array([1/4.,0.,0.])
    w = np.linspace(0, 24, 241)
    
    df = DF(calc='Al.gpw',
            q=q,
            w=w,
            eta=0.2,
            ecut=50,
            hilbert_trans=False)
    df.get_EELS_spectrum(filename='Al_df.dat')
    df.write('Al.pckl')
    df.check_sum_rule()


if check_spectrum:
    d = np.loadtxt('Al_bse.dat')[:,2] 
    wpeak = 16.4 
    Nw = 164
    if d[Nw] > d[Nw-1] and d[Nw] > d[Nw+1]:
        pass
    else:
        raise ValueError('Plasmon peak not correct ! ')
    
    if np.abs(d[Nw] - 27.4958893542) > 1e-5:
        print(d[Nw])
        raise ValueError('Please check spectrum strength ! ')

    d2 = np.loadtxt('Al_df.dat')
    if np.abs(d[:240] - d2[:240, 2]).sum() > 0.003:
        raise ValueError('Please compare two spectrum')
