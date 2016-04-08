from gpaw.tddft import *
from ase import Atoms
from gpaw import GPAW
from sys import argv
from gpaw.tddft.lcao_tddft import LCAOTDDFT
import gpaw.lrtddft as lrtddft
from sys import exit
from gpaw.mpi import world
xc = 'oldLDA'
c = +1
h = 0.4
b = 'dzp'
sy = 'Na2'
positions = [[0.00,0.00,0.00],[0.00,0.00,2.00]]
atoms = Atoms(symbols=sy, positions = positions)
atoms.center(vacuum=3)

# LCAO-RT-TDDFT
calc = GPAW(mode='lcao', nbands=1, xc=xc, h=h, basis=b, dtype=complex, charge=c, width=0, convergence={'density':1e-8})
atoms.set_calculator(calc)
atoms.get_potential_energy()
calc.write('Na2.gpw','all')
del calc

calc = LCAOTDDFT('Na2.gpw')
dmfile = sy+'_lcao_restart_'+b+'_rt_z.dm'+str(world.size)
specfile = sy+'_lcao_restart_'+b+'_rt_z.spectrum'+str(world.size)
calc.absorption_kick([0.0,0,0.001])
calc.propagate(10, 20, dmfile)
if world.rank == 0:
    photoabsorption_spectrum(dmfile, specfile)

