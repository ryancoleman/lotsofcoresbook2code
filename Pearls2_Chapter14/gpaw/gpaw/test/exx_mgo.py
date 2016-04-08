from __future__ import print_function
import numpy as np
from ase.units import Ha
from ase.dft.kpoints import monkhorst_pack
from ase.parallel import paropen
from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from gpaw.wavefunctions.pw import PW
from gpaw.mpi import size
from gpaw.xc.hybridg import HybridXC
from ase.io import read

for nk in (4,6,8):

    kpts = monkhorst_pack((nk,nk,nk))
    kshift = 1./(2*nk)
    kpts += np.array([kshift, kshift, kshift])

    atoms = bulk('MgO', 'rocksalt', a = 4.212)

    calc = GPAW(mode=PW(600),dtype=complex, basis='dzp', xc='PBE', maxiter=300,
            txt='gs_%s.txt'%(nk), kpts=kpts,parallel={'band':1},
                        occupations=FermiDirac(0.01))
    
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    E = calc.get_potential_energy()
    E_exx = E + calc.get_xc_difference(HybridXC('EXX', method='acdf'))

    print(nk, E, E_exx)


 
