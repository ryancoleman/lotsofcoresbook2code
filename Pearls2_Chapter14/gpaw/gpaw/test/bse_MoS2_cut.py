import os
import numpy as np
from ase import *
from ase.lattice.hexagonal import Hexagonal
from gpaw import *
from gpaw.response.bse import BSE
from gpaw.mpi import rank

calc = GPAW(h=0.2,
            xc='GLLBSC',
            nbands=20,
            setups={'Mo': '6'},
            occupations=FermiDirac(0.001),
            convergence={'bands': -5},
            kpts=(9,9,1))

a = 3.1604
c = 10.0

cell = Hexagonal(symbol='Mo', latticeconstant={'a':a, 'c':c}).get_cell()
layer = Atoms(symbols='MoS2', cell=cell, pbc=(1,1,1),
              scaled_positions=[(0, 0, 0),
                                (2/3., 1/3., 0.3),
                                (2/3., 1/3., -0.3)])

pos = layer.get_positions()
pos[1][2] = pos[0][2] + 3.172/2
pos[2][2] = pos[0][2] - 3.172/2
layer.set_positions(pos)
layer.set_calculator(calc)
layer.get_potential_energy()
response = calc.hamiltonian.xc.xcs['RESPONSE']
response.calculate_delta_xc()
E_ks, dis = response.calculate_delta_xc_perturbation()
    

bse = BSE(calc,
          w=np.linspace(0., 5., 501),
          q=np.array([0.0001, 0., 0.]),
          optical_limit=True,
          ecut=10,
          eta=0.02,
          nv=np.array([8,9]),
          nc=np.array([9,10]),
          eshift=dis,
          nbands=15,
          mode='BSE',
          vcut='2D', 
          )

bse.get_dielectric_function('bse_cut.dat')

if rank == 0 and os.path.isfile('phi_qaGp'):
    os.remove('phi_qaGp')
    
    d = np.loadtxt('bse_cut.dat')
    Nw = 88
    if d[Nw,2] > d[Nw-1,2] and d[Nw,2] > d[Nw+1,2]:
        pass
    else:
        raise ValueError('Absorption peak not correct ! ')
