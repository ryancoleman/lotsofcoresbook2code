#!/usr/bin/env python
#PBS -N NiC8
#PBS -e NiC8.err
#PBS -o NiC8.log
#PBS -m ae
#PBS -q verylong
#PBS -l nodes=4:ppn=8:xeon5570

from gpaw import *
from ase import *
from gpaw.transport.calculator import Transport
from ase.lattice.surface import fcc100
from gpaw.transport.tools import sort_atoms

system = read('NiC8.traj')
system.set_pbc([1,1,1])
sort_atoms(system)
system.center()
lead_cell = np.diag(system.cell)
lead_cell[2] = 7.040

pl_atoms1 = range(36)     
pl_atoms2 = range(98,134) 
pl_cell1 = lead_cell
pl_cell2 = pl_cell1      

magmoms = np.zeros([134])
magmoms[:54]=0.7
magmoms[-54:] = -0.7
system.set_initial_magnetic_moments(magmoms)

t = Transport( h=0.2,     
               xc='RPBE',
               basis={'Ni': 'szp', 'H': 'szp', 'C': 'szp', 'S':'szp'},
               kpts=(2, 2, 1),
               occupations=FermiDirac(0.2),
               mode='lcao',
               txt='NiC8.txt',
               buffer_guess=True,
               lead_guess=True,
               spinpol=True,
               guess_steps=80,
               beta_guess=0.003,
               alpha=0.1,
               poissonsolver=PoissonSolver(nn=2),
               mixer=MixerSum(0.005, 5, weight=100.0),
               extra_density=True,
               pl_atoms=[pl_atoms1, pl_atoms2],
               pl_cells=[pl_cell1, pl_cell2],
               pl_kpts=[2 , 2 , 15],
               edge_atoms=[[ 0, 35],[0 , 133]],
               mol_atoms=range(36, 98),
               nleadlayers=[1,1]) 
system.set_calculator(t)
t.calculate_iv()

t = Transport( h=0.2,     
               xc='RPBE',
               basis={'Ni': 'szp', 'H': 'szp', 'C': 'szp', 'S':'szp'},
               kpts=(12, 12, 1),
               occupations=FermiDirac(0.2),
               parallel={'domain':(1,1,1)},
               mode='lcao',
               txt='NiC8.txt',
               buffer_guess=True,
               lead_guess=True,
               spinpol=True,
               guess_steps=80,
               beta_guess=0.003,
               alpha=0.1,
               poissonsolver=PoissonSolver(nn=2),
               mixer=MixerSum(0.005, 5, weight=100.0),
               extra_density=True,
               analysis_mode=True,
               pl_atoms=[pl_atoms1, pl_atoms2],
               pl_cells=[pl_cell1, pl_cell2],
               pl_kpts=[12 , 12 , 15],
               edge_atoms=[[ 0, 35],[0 , 133]],
               mol_atoms=range(36, 98),
               nleadlayers=[1,1]) 
system.set_calculator(t)
t.analysis(1)

