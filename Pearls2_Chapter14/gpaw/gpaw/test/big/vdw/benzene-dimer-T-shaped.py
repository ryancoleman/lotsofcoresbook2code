from __future__ import print_function
import numpy as np
from ase.structure import molecule
from ase.constraints import FixedPlane
from ase.optimize import QuasiNewton
from gpaw import GPAW, FermiDirac
from gpaw.test import equal

#----------------------------------
# Initialization
molname = 'benzene-mol'
dimername = 'benzene-dimer'
f = open('benzene-dimer-T-shape.dat', 'w')
h = 0.18
xc = 'vdW-DF'

#-------------------------------------
# relaxation of the benzene molecule
benz = molecule('C6H6')
benz.set_pbc(False)
tags = np.zeros_like(benz)
benz.set_tags(tags)
benz.center(vacuum=4.0)
cell = benz.get_cell()

calc = GPAW(nbands=-1,
            h=h,
            xc=xc,
            occupations=FermiDirac(0.0),
            txt=molname+'_relax.txt')
benz.set_calculator(calc)

# qn constraint
for i in range(len(benz)):
    plane = FixedPlane(i, (0, 0, 1))
    benz.set_constraint(plane)

qn = QuasiNewton(benz,
                 logfile=molname + '_relax.log',
                 trajectory=molname + '_relax.traj')
qn.run(fmax=0.01)

e_mol = benz.get_potential_energy()
del calc

#-------------------------------------
# mapping out the benzene dimer (T-shaped) intermolecular distance

e_array = np.zeros(20)
d_array = np.zeros(20)

k = 0
for i in np.linspace(-6, 6, 20):
    k_str = str(k)
    z = 6.0 + i * 0.3
    BENZ = benz.copy()
    dimer = BENZ.copy()
    tags = np.ones_like(dimer)
    dimer.set_tags(tags)
    BENZ.rotate('x', np.pi / 2, center='COM')
    BENZ.translate([0, 0, z])
    dimer.extend(BENZ)
    dimer.set_cell([cell[0, 0], cell[1, 1], cell[2, 2] + 8])
    dimer.center()
    dimer.set_pbc(False)
    pos = dimer.get_positions()
    d = pos[21, 2] - pos[0, 2]
    
    calc = GPAW(nbands=-2,
                h=h,
                xc=xc,
                occupations=FermiDirac(0.0),
                txt=dimername + '_' + k_str + '.txt')
    dimer.set_calculator(calc)
    e_dimer = dimer.get_potential_energy()
    del calc
    
    # interaction energy
    e_int = e_dimer - 2 * e_mol
    e_array[k] = e_int
    d_array[k] = d
    
    print(str(round(d, 3)), e_int, file=f)
    f.flush()
    k += 1

# E_int-curve minimum
e_0 = 100.
d_0 = 0.
for i in range(len(e_array)):
    if e_array[i] < e_0:
        e_0 = e_array[i]
        d_0 = d_array[i]
print('****************', file=f)
print('Minimum (E_int,d):', e_0, d_0, file=f) 
f.close()

equal(e_0 , -0.11, 0.01)
equal(d_0 , 2.86, 0.05)
