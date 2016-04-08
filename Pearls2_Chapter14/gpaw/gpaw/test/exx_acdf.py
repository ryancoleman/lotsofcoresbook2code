from __future__ import print_function
from ase import *
from ase.dft.kpoints import monkhorst_pack
from ase.lattice import bulk
from gpaw import *
from gpaw.test import equal
import numpy as np
from gpaw.mpi import serial_comm, size, rank, world

atoms = bulk('Al', 'fcc')
kpts = monkhorst_pack((4,4,4))
kpts += np.array([1/8., 1/8., 1/8.])

calc = GPAW(h=0.18,
            kpts=kpts,
            occupations=FermiDirac(0.1))
atoms.set_calculator(calc)
E = atoms.get_potential_energy()
calc.write('Al.gpw','all') 

calc = GPAW('Al.gpw',txt=None)
E = calc.get_potential_energy()

from gpaw.xc.hybridk import HybridXC
exx = HybridXC('EXX',acdf=True)
E_k = E + calc.get_xc_difference(exx)

if 0:#size == 1:
    calc = GPAW('Al.gpw',txt=None, communicator=serial_comm)
    from gpaw.xc.hybridq import HybridXC
    exx = HybridXC('EXX')
    E_q = E + calc.get_xc_difference(exx)

    print(E_q, E_k)
    equal(E_q, E_k, 0.001)
equal(E_k, -15.6, 0.1)
