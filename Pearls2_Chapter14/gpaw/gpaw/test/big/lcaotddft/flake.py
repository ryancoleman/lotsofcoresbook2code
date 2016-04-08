# Standard magic
from ase import Atoms
from ase.io import read
from ase.parallel import parprint, paropen
from gpaw import GPAW, restart, Mixer, PoissonSolver, FermiDirac
from gpaw.tddft import *
from gpaw.tddft.lcao_tddft import LCAOTDDFT

basis = 'dzp'
spacing = 0.32
title = 'flake'
vacuum = spacing*8*4/2
time_step = 5.0
kick = [0.001, 0.000, 0.000 ]

from os import path
if path.exists('flake.gpw'):
    calc = LCAOTDDFT('flake.gpw')
    atoms = calc.get_atoms()
else:
    calc = LCAOTDDFT(mode='lcao', h=spacing, basis=basis, 
                        nbands=423, width=0, 
                        mixer=Mixer(0.05, 5, weight=100.0),
                       poissonsolver=PoissonSolver(eps=1e-12))
    atoms = read('%s.xyz' % (title));
    atoms.set_pbc((False, False, False))
    atoms.center(vacuum=vacuum)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    gs_calc.write('flake.gpw', 'all')
  
maxiterations = 24000/time_step
fname0 = 'flake_dm.dat'
fname2 = 'flake_spectrum.dat'
  
calc.absorption_kick(kick)
calc.propagate(time_step, maxiterations, fname0)
 
#photoabsorption_spectrum(fname0, fname2, e_min=0.0, e_max=40.0, delta_e=0.02, width=0.05)
