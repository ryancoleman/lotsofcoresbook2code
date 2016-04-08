import numpy as np
from ase.units import Ha
from ase.dft.kpoints import monkhorst_pack
from ase.parallel import paropen
from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from gpaw.wavefunctions.pw import PW
from gpaw.mpi import size

kpts = monkhorst_pack((10,10,10))
kpts += np.array([1/20., 1/20., 1/20.])

bulk = bulk('Na', 'bcc', a=4.23)

tag = 'Nabulk'

if 1:
    ecut = 350
    calc = GPAW(mode=PW(ecut),dtype=complex, basis='dzp', kpts=kpts, xc='PBE',
                eigensolver='rmm-diis',
                parallel={'band': size}, txt='gs_occ_%s.txt' %(tag), nbands=4,
                occupations=FermiDirac(0.01), setups={'Na': '1'},
                )
    bulk.set_calculator(calc)
    bulk.get_potential_energy()
    calc.write('gs_occ_%s.gpw' %(tag))

if 1:
    calc = GPAW('gs_occ_%s.gpw' %(tag),txt='gs_%s.txt'%(tag), parallel={'band': 1, 'domain':1})
    calc.diagonalize_full_hamiltonian(nbands=520)
    calc.write('gs_%s.gpw' %(tag), 'all')

