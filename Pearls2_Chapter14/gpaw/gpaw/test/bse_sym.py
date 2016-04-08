from __future__ import print_function
import os
import numpy as np
from ase.lattice import bulk
from gpaw import GPAW
from gpaw.response.df0 import DF
from ase.dft.kpoints import monkhorst_pack
from gpaw.response.bse import BSE
from gpaw.mpi import rank, size

# generate kmesh
kpts =(2,2,2)
bzk_kc = monkhorst_pack(kpts)

shift_c = []
for Nk in kpts:
    if Nk % 2 == 0:
        shift_c.append(0.5 / Nk)
    else:
        shift_c.append(0.)

atoms = bulk('Si', 'diamond', a=5.431) 

kpts1 = bzk_kc # not Gamma centered
kpts2 = bzk_kc + shift_c # Gamma centered

for kpts in (kpts2,):

    calc = GPAW(h=0.20, kpts=kpts)      
    atoms.set_calculator(calc)               
    atoms.get_potential_energy()
    calc.write('Si.gpw','all')
    
    # no symmetry BSE
    eshift = 0.8
    bse = BSE('Si.gpw',
              w=np.linspace(0,10,201),
              q=np.array([0.0001, 0, 0.0]),
              optical_limit=True,
              ecut=150.,
              nc=np.array([4,6]),
              nv=np.array([2,4]),
              eshift=eshift,
              nbands=8,
              qsymm=False)
    bse.get_dielectric_function('bse_nosymm.dat')

    if rank == 0 and os.path.isfile('phi_qaGp'):
        os.remove('phi_qaGp')

    
    # with symmetry BSE
    eshift = 0.8
    bse = BSE('Si.gpw',
              w=np.linspace(0,10,201),
              q=np.array([0.0001,0,0.0]),
              optical_limit=True,
              ecut=150.,
              nc=np.array([4,6]),
              nv=np.array([2,4]),
              eshift=eshift,
              nbands=8,
              qsymm=True)
    bse.get_dielectric_function('bse_symm.dat')

    if rank == 0 and os.path.isfile('phi_qaGp'):
        os.remove('phi_qaGp')

check = 1
if check:
    d1 = np.loadtxt('bse_nosymm.dat')
    d2 = np.loadtxt('bse_symm.dat')
    print(np.abs(d1[:,2] - d2[:,2]).max(), np.abs(d1[:,2] - d2[:,2]).sum())
    assert np.abs(np.abs(d1[:,2] - d2[:,2]).max() - 0.014775742) < 1e-2
    assert np.abs(np.abs(d1[:,2] - d2[:,2]).sum() - 0.210880672212) < 1e-1
