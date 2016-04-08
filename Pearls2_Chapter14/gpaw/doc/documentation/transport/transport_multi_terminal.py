from ase import Atoms, view
from gpaw import FermiDirac, Mixer
from gpaw.transport.calculator import Transport 
from gpaw.atom.basis import BasisMaker
from gpaw.poisson import PoissonSolver
import pickle
import numpy as np

Ly=7.0
Lx= 30.
Lz= 30.
C_C=1.42
C_H=1.1
basis = 'sz'
atoms = Atoms('C18H3', pbc=(0, 0, 0), cell=[Lx, Ly, Lz])
atoms.positions[:, 1] = 0
for i in range(6):
    atoms.positions[i, 0] = np.cos(np.pi/3*i)*C_C
    atoms.positions[i, 2] = np.sin(np.pi/3*i)*C_C
for i in range(6,10):
    atoms.positions[i, 0] = (i + 2 - 6) * atoms.positions[0, 0] 
    atoms.positions[i, 2] = (i + 2 - 6) * atoms.positions[0, 2]
for i in range(10,14):
    atoms.positions[i, 0] = (i + 2 - 10) * atoms.positions[2, 0] 
    atoms.positions[i, 2] = (i + 2 - 10) * atoms.positions[2, 2]
for i in range(14,18):
    atoms.positions[i, 0] = (i + 2 - 14) * atoms.positions[4, 0] 
    atoms.positions[i, 2] = (i + 2 - 14) * atoms.positions[4, 2]
for i in range(18,21):
    atoms.positions[i, 0] = atoms.positions[(i - 18)*2 + 1, 0] * (C_H + C_C) / C_C
    atoms.positions[i, 2] = atoms.positions[(i - 18)*2 + 1, 2] * (C_H + C_C) / C_C

#sort_atoms(atoms)
atoms.rotate('x', 'z')
atoms.rotate('y', np.pi)
atoms.center()

eatoms = Atoms('C18H3C12', pbc=(0, 0, 0), cell=[Lx, Ly, Lz])
eatoms.positions[:, 1] = 0
for i in range(6):
    eatoms.positions[i, 0] = np.cos(np.pi/3*i)*C_C
    eatoms.positions[i, 2] = np.sin(np.pi/3*i)*C_C
for i in range(6,10):
    eatoms.positions[i, 0] = (i + 2 - 6) * eatoms.positions[0, 0] 
    eatoms.positions[i, 2] = (i + 2 - 6) * eatoms.positions[0, 2]
for i in range(10,14):
    eatoms.positions[i, 0] = (i + 2 - 10) * eatoms.positions[2, 0] 
    eatoms.positions[i, 2] = (i + 2 - 10) * eatoms.positions[2, 2]
for i in range(14,18):
    eatoms.positions[i, 0] = (i + 2 - 14) * eatoms.positions[4, 0] 
    eatoms.positions[i, 2] = (i + 2 - 14) * eatoms.positions[4, 2]
for i in range(18,21):
    eatoms.positions[i, 0] = eatoms.positions[(i - 18)*2 + 1, 0] * (C_H + C_C) / C_C
    eatoms.positions[i, 2] = eatoms.positions[(i - 18)*2 + 1, 2] * (C_H + C_C) / C_C

for i in range(21,25):
    eatoms.positions[i, 0] = (i + 2 - 17) * eatoms.positions[0, 0] 
    eatoms.positions[i, 2] = (i + 2 - 17) * eatoms.positions[0, 2]
for i in range(25,29):
    eatoms.positions[i, 0] = (i + 2 - 21) * eatoms.positions[2, 0] 
    eatoms.positions[i, 2] = (i + 2 - 21) * eatoms.positions[2, 2]
for i in range(29,33):
    eatoms.positions[i, 0] = (i + 2 - 25) * eatoms.positions[4, 0] 
    eatoms.positions[i, 2] = (i + 2 - 25) * eatoms.positions[4, 2]

#sort_atoms(atoms)
eatoms.rotate('x', 'z')
eatoms.rotate('y', np.pi)
eatoms.center()

pl_atoms1 = range(6,10)     
pl_atoms2 = range(10,14)
pl_atoms3 = range(14,18)

pl_cell1 = [Ly, Ly, 4*C_C]
pl_cell2 = pl_cell1
pl_cell3 = pl_cell1

t = Transport(h=0.25,
              xc='LDA',
              basis='sz',
              occupations=FermiDirac(0.2),
              kpts=(1,1,1),
              mode='lcao',
              poissonsolver=PoissonSolver(nn=2, relax='GS'),
              txt='Na_lcao.txt',
              mixer=Mixer(0.1, 5, weight=100.0),
              guess_steps=10,
              LR_leads=False,
              bias=[0]*3,
              fixed_boundary=False,
              identical_leads=True,
              pl_atoms=[pl_atoms1, pl_atoms2, pl_atoms3],
              pl_cells=[pl_cell1, pl_cell2, pl_cell3],
              pl_kpts=(1,1,15),
              save_file=True,
              use_buffer=True,
              neutral=False,
              normalize_density=False,
              extended_atoms=eatoms,
              multi_lead_directions=[[6, 9], [10, 13], [14,17]],
              nleadlayers=[1,1,1],
              edge_atoms=[[0, 0, 0], [6, 10, 14]],
              mol_atoms=range(6)+ range(18,21))
atoms.set_calculator(t)

lp = []
for i in range(3):
    for j in range(i+1, 3):
        lp.append([i,j])

t.set_analysis_parameters(lead_pairs=lp)
t.negf_prepare()
t.get_selfconsistent_hamiltonian()

