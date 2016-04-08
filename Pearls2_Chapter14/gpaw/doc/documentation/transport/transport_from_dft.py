from ase import Atoms
from gpaw.transport.calculator import Transport 
from gpaw.atom.basis import BasisMaker
from gpaw.occupations import FermiDirac
from gpaw.poisson import PoissonSolver
from gpaw.mixer import Mixer
from gpaw import GPAW
import pickle

a = 3.6
L = 7.00

basis = BasisMaker('Na').generate(1, 1, energysplit=0.3)
#Lead1 and Lead2 can be different
Lead1 = Atoms('Na4', pbc=(1, 1, 1), cell=[L, L, 4 * a])
Lead1.positions[:4, 2] = [i * a for i in range(4)]
Lead1.positions[:, :2] = L / 2.
Lead1.center()

Lead1.set_calculator(GPAW(h=0.3,
                          xc='LDA',
                          basis={'Na': basis},
                          kpts=(1,1,5),
                          occupations=FermiDirac(0.1),
                          mode='lcao',
                          mixer=Mixer(0.1, 5, weight=100.0),
                          txt='Lead1.txt'))
Lead1.get_potential_energy()

Lead2 = Lead1.copy()
Lead2.set_calculator(GPAW(h=0.3,
                          xc='LDA',
                          basis={'Na': basis},
                          kpts=(1,1,5),
                          occupations=FermiDirac(0.1),
                          mode='lcao',
                          mixer=Mixer(0.1, 5, weight=100.0),
                          txt='Lead2.txt'))
Lead2.get_potential_energy()

Device = Atoms('Na9', pbc=(1, 1, 1), cell=[L, L, 9 * a])
Device.positions[:9, 2] = [i * a for i in range(9)]
Device.positions[:, :2] = L / 2.
Device.center()

Device.set_calculator(GPAW(h=0.3,
                           xc='LDA',
                           basis={'Na': basis},
                           kpts=(1,1,1),
                           occupations=FermiDirac(0.1),
                           mode='lcao',
                           mixer=Mixer(0.1, 5, weight=100.0),
                           txt='Device.txt'))
Device.get_potential_energy()

from gpaw.transport.tools import save_bias_data_file
save_bias_data_file(Lead1, Lead2, Device)

system = Device.copy()
pl_cell1 = (L, L, 4 * a) 
pl_cell2 = pl_cell1
t = Transport(h=0.3,
              xc='LDA',
              basis={'Na': basis},
              kpts=(1,1,1),
              occupations=FermiDirac(0.1),
              mode='lcao',
              txt='Na_lcao.txt',
              mixer=Mixer(0.1, 5, weight=100.0),
              pl_atoms=[range(4), range(5,9)],
              pl_cells=[pl_cell1, pl_cell2],
              leads=[Lead1,Lead2],
              pl_kpts=(1,1,5),
              analysis_mode=True)
system.set_calculator(t)
t.analysis(1)
