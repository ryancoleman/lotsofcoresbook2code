from __future__ import print_function
import os
import sys

from ase import Atom, Atoms
from ase.units import Bohr, Hartree
from ase.io.cube import write_cube
from gpaw import GPAW
from gpaw.test import equal

from gpaw.cluster import Cluster
from gpaw.point_charges import PointCharges
from gpaw.external_potential import ConstantPotential
from gpaw.mpi import world

cp = ConstantPotential()

sc = 2.9
R=0.7 # approx. experimental bond length
R=1.
a = 2 * sc
c = 3 * sc
at='H'
H2 = Atoms([Atom(at, (a/2, a/2, (c-R)/2)),
            Atom(at, (a/2, a/2, (c+R)/2))],
           cell=(a,a,c), pbc=False)
print(at, 'dimer')
nelectrons = 2 * H2[0].number

txt = None
#txt = '-'

# load point charges
fname = 'pc.xyz'
if world.rank == 0:
    f = open('i' + fname, 'w')
    print("""1

X 0 0 100 -0.5""", file=f)
    f.close()
world.barrier()

ex = PointCharges()
ex.read('i' + fname)
ex.write('o' + fname)

convergence = {'eigenstates':1.e-4*40*1.5**3, 'density':1.e-2, 'energy':0.1}

# without potential
if True:
    if txt:
        print('\n################## no potential')
    c00 = GPAW(h=0.3, nbands=-1,
               convergence=convergence,
               txt=txt)
    c00.calculate(H2)
    eps00_n = c00.get_eigenvalues()

# 0 potential
if True:
    if txt:
        print('\n################## 0 potential')
    cp0 = ConstantPotential(0.0)
    c01 = GPAW(h=0.3, nbands=-2, external=cp0, 
               convergence=convergence,
               txt=txt)
    c01.calculate(H2)

# 1 potential
if True:
    if txt:
        print('################## 1 potential')
    cp1 = ConstantPotential(-1.0/Hartree)
    c1 = GPAW(h=0.3, nbands=-2, external=cp1, 
              convergence=convergence,
              txt=txt)
    c1.calculate(H2)

for i in range(c00.get_number_of_bands()):
    f00 = c00.get_occupation_numbers()[i]
    if f00 > 0.01:
        e00 = c00.get_eigenvalues()[i]
        e1 = c1.get_eigenvalues()[i]
        print('Eigenvalues no pot, expected, error=', e00, e1 + 1, e00 - e1 - 1)
        equal(e00, e1 + 1., 0.007)

E_c00 = c00.get_potential_energy()
niter_c00 = c00.get_number_of_iterations()

E_c1 = c1.get_potential_energy()
niter_c1 = c1.get_number_of_iterations()

DeltaE = E_c00 - E_c1
print('Energy diff, expected, error=', DeltaE, nelectrons, DeltaE - nelectrons)
equal(DeltaE, nelectrons, 0.002)

energy_tolerance = 0.00001
niter_tolerance = 0
#equal(E_c00, 10.4409370467, energy_tolerance) # svnversion 5252
#equal(niter_c00, 14, niter_tolerance) # svnversion 5252
#equal(E_c1, -11.5590572387, energy_tolerance) # svnversion 5252
#equal(niter_c1, 14, niter_tolerance) # svnversion 5252
