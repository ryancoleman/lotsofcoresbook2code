import numpy as np
from ase.io import read
from gpaw import FermiDirac, PoissonSolver, Mixer
from gpaw.transport.calculator import Transport

system = read('BDT.traj', -1)
pl_atoms1 = range(36)     
pl_atoms2 = range(72, 108) 
pl_cell1 = np.diag(system.cell)
pl_cell1[2] = 7.067
pl_cell2 = pl_cell1      

t = Transport(h=0.2,     
              xc='PBE',
              basis={'Au': 'sz(dzp)', 'H': 'dzp', 'C': 'dzp', 'S':'dzp'},
              kpts=(1, 1, 1),
              occupations=FermiDirac(0.1),
              mode='lcao',
              txt='ABA.txt',
              poissonsolver=PoissonSolver(nn=2),
              mixer=Mixer(0.1, 5, weight=100.0),
              pl_atoms=[pl_atoms1, pl_atoms2],
              pl_cells=[pl_cell1, pl_cell2],
              pl_kpts=[1, 1, 15],
              extra_density=True,
              analysis_data_list=['tc'],
              edge_atoms=[[ 0, 35],[0 , 107]],
              mol_atoms=range(36, 72))
system.set_calculator(t)
t.calculate_iv(3.0, 2)
