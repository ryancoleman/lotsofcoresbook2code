from __future__ import print_function
from ase import *
from ase.dft.kpoints import monkhorst_pack
from ase.lattice import bulk
from gpaw import *
from gpaw.test import equal
import numpy as np

a0 = 5.43
cell = bulk('Si', 'fcc', a=a0).get_cell()
Si = Atoms('Si2', cell=cell, pbc=True,
           scaled_positions=((0,0,0), (0.25,0.25,0.25)))

kpts = monkhorst_pack((2,2,2))
kpts += np.array([1/4., 1/4., 1/4.])

calc = GPAW(h=0.18,
            kpts=kpts,
#            parallel={'domain':1}, idiotproof=False,
            occupations=FermiDirac(0.001))
Si.set_calculator(calc)
E = Si.get_potential_energy()

#from gpaw.xc.hybridq import HybridXC
#exx = HybridXC('EXX')
#E_q = E + calc.get_xc_difference(exx)

from gpaw.xc.hybridk import HybridXC
exx = HybridXC('EXX', acdf=True, etotflag=True)
E_k1 = E + calc.get_xc_difference(exx)

from gpaw.xc.hybridk import HybridXC
exx = HybridXC('EXX', acdf=False, etotflag=True)
E_k2 = E + calc.get_xc_difference(exx)

#print 'Hartree-Fock ACDF method (hybridq.py)   :', E_q
print('Hartree-Fock ACDF method (hybridk.py)   :', E_k1)
print('Hartree-Fock Standard method            :', E_k2)

equal(E_k2, E_k1, 0.001)
#equal(E_q, E_k1, 0.001)
equal(E_k1, -27.71, 0.01)
