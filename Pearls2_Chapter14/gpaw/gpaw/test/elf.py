from __future__ import print_function
import sys
import numpy as np
from ase.structure import molecule
from ase.io import write
from ase.parallel import parprint
from gpaw import GPAW, restart
from gpaw.elf import ELF
from gpaw.test import equal
from gpaw.mpi import rank, world

atoms = molecule('CO')
atoms.center(2.0)

txt=sys.stdout
txt=None

try:
    atoms, calc = restart('CO.gpw', txt=txt)
    energy = atoms.get_potential_energy()
except:
    calc = GPAW(h=0.24, txt=txt)
    atoms.set_calculator(calc)
    energy = atoms.get_potential_energy()
    calc.write('CO.gpw', 'all')

elf = ELF(calc)
elf.update()
elf_G = elf.get_electronic_localization_function(gridrefinement=1)
elf_g = elf.get_electronic_localization_function(gridrefinement=2)

nt_G = calc.density.nt_sG[0]
taut_G = elf.taut_sG[0]
nt_grad2_G = elf.nt_grad2_sG[0]
nt_grad2_g = elf.nt_grad2_sg[0]

# integrate the CO bond
if rank == 0:
    # bond area
    x0 = (atoms.positions[0][0] - 1.0)/atoms.get_cell()[0,0]
    x1 = 1 - x0
    y0 = (atoms.positions[0][1] -1.0)/atoms.get_cell()[1,1]
    y1 = 1 - y0
    z0 = atoms.positions[1][2]/atoms.get_cell()[2,2]
    z1 = atoms.positions[0][2]/atoms.get_cell()[2,2]
    gd = calc.wfs.gd
    Gx0, Gx1 = int(gd.N_c[0]*x0), int(gd.N_c[0] * x1)
    Gy0, Gy1 = int(gd.N_c[1]*y0), int(gd.N_c[1] * y1)
    Gz0, Gz1 = int(gd.N_c[2]*z0), int(gd.N_c[2] * z1)
    finegd = calc.density.finegd
    gx0, gx1 = int(finegd.N_c[0]*x0), int(finegd.N_c[0] * x1)
    gy0, gy1 = int(finegd.N_c[1]*y0), int(finegd.N_c[1] * y1)
    gz0, gz1 = int(finegd.N_c[2]*z0), int(finegd.N_c[2] * z1)
    int1 = elf_G[Gx0:Gx1,Gy0:Gy1,Gz0:Gz1].sum() * gd.dv
    int2 = elf_g[gx0:gx1,gy0:gy1,gz0:gz1].sum() * finegd.dv
    parprint("Ints", int1, int2)
    parprint("Min, max G", np.min(elf_G), np.max(elf_G))
    parprint("Min, max g", np.min(elf_g), np.max(elf_g))
#   The tested values (< r7887) do not seem to be correct  
#    equal(int1, 14.8078, 0.0001)
#    equal(int2, 13.0331, 0.0001)

# check spin-polarized version
try:
    atoms, calc = restart('COspin.gpw', txt=None,
                          parallel={'domain': world.size})
    energy_spinpol = atoms.get_potential_energy()
except:
    calc.set(spinpol=True,
             parallel={'domain': world.size})
    energy_spinpol = atoms.get_potential_energy()
    calc.write('COspin.gpw', 'all')

def check_diff(g1, g2, gd, txt):
#    print rank, txt, "shapes", g1.shape, g2.shape
    intd = gd.integrate(np.abs(g1 - g2))
    parprint(txt, 'integrated diff=', intd, end='')
    maxd = np.max(np.abs(g1 - g2))
    parprint('max diff=', maxd)
    equal(intd, 0, 1.e-8) 
    equal(maxd, 0, 1.e-9) 

nt_spinpol_G = calc.density.nt_sG.sum(axis=0)
check_diff(nt_G, nt_spinpol_G, elf.finegd, 'nt_G')
   
equal(energy, energy_spinpol, 0.0001)

elf_spinpol = ELF(calc)
elf_spinpol.update()
elf_spinpol_G = elf_spinpol.get_electronic_localization_function(gridrefinement=1)
elf_spinpol_g = elf_spinpol.get_electronic_localization_function(gridrefinement=2)
taut_spinpol_G = elf_spinpol.taut_sG.sum(axis=0)
check_diff(taut_G, taut_spinpol_G, elf.gd, 'taut_G')

nt_grad2_spinpol_G = 2 * elf_spinpol.nt_grad2_sG.sum(axis=0)
check_diff(nt_grad2_G, nt_grad2_spinpol_G, elf.gd, 'nt_grad2_G')

nt_grad2_spinpol_g = 2 * elf_spinpol.nt_grad2_sg.sum(axis=0)
check_diff(nt_grad2_g, nt_grad2_spinpol_g, elf.finegd, 'nt_grad2_g')

check_diff(elf_G, elf_spinpol_G, elf.gd, 'elf_G')
check_diff(elf_g, elf_spinpol_g, elf.finegd, 'elf_g')
