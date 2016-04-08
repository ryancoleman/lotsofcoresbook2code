from gpaw import *
from gpaw.transport.calculator import Transport 
from gpaw.atom.basis import BasisMaker
from gpaw.poisson import PoissonSolver
from ase.optimize import QuasiNewton
from ase import Atoms

a = 3.6
L = 7.00

basis = 'sz'

atoms = Atoms('Na12', pbc=(0, 0, 1), cell=[L, L, 12 * a + 2.0])
atoms.positions[:6, 2] = [i * a for i in range(6)]
atoms.positions[6:12, 2] = [i * a + 2.0 for i in range(6, 12)]

atoms.positions[:, :2] = L / 2.
atoms.center()

pl_atoms1 = range(4)     
pl_atoms2 = range(8, 12)
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
              guess_steps=10,
              fix_contour=True,
              pl_atoms=[pl_atoms1, pl_atoms2],
              pl_cells=[pl_cell1, pl_cell2],
              pl_kpts=(1,1,16),
              fixed_boundary=False,
              analysis_data_list=['tc'],
              edge_atoms=[[0, 3], [0, 11]],
              mol_atoms=range(4, 8))
atoms.set_calculator(t)
t.calculate_to_bias(0.0, 1)
qn=QuasiNewton(atoms, trajectory='transport.traj')
qn.run(fmax=0.05)

