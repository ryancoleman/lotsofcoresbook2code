from __future__ import print_function
import pickle
from ase import Atoms
from gpaw.transport.calculator import Transport 
from gpaw.atom.basis import BasisMaker
from gpaw.poisson import PoissonSolver
from gpaw.mixer import Mixer
from gpaw.occupations import FermiDirac

a = 3.6
L = 7.00

basis = BasisMaker('Na').generate(1, 1, energysplit=0.3)

atoms = Atoms('Na9',
              cell=(L, L, 9 * a),
              pbc=(0, 0, 1))
              
atoms.positions[:9, 2] = [i * a for i in range(9)]
atoms.positions[:, :2] = L / 2.
atoms.center()

pl_atoms1 = range(4)     
pl_atoms2 = range(5, 9)
pl_cell1 = (L, L, 4 * a) 
pl_cell2 = pl_cell1

t = Transport(h=0.3,
              xc='LDA',
              basis={'Na': basis},
              kpts=(1,1,1),
              occupations=FermiDirac(0.1),
              mode='lcao',
              poissonsolver=PoissonSolver(nn=2, relax='GS'),
              txt='Na_lcao.txt',
              mixer=Mixer(0.1, 5, weight=100.0),
              pl_atoms=[pl_atoms1, pl_atoms2],
              pl_cells=[pl_cell1, pl_cell2],
              pl_kpts=(1,1,5),
              edge_atoms=[[0, 3], [0, 8]],
              mol_atoms=[4],
              guess_steps=1,
              fixed_boundary=False,
              foot_print=False)

atoms.set_calculator(t)
t.calculate_iv(0.5, 2)

