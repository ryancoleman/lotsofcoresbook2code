from __future__ import print_function
import numpy as np
from ase.lattice import bulk
from ase.dft.kpoints import monkhorst_pack
from ase.parallel import paropen

from gpaw import GPAW, PW
from gpaw.mpi import size, rank, world, serial_comm
from gpaw.xc.tools import vxc
from gpaw.xc.hybridg import HybridXC


mgo = bulk('MgO', 'rocksalt', a=4.189)
if rank < 3:
    comm = world.new_communicator(np.arange(min(3, size)))
else:
    comm = world.new_communicator(np.array((rank,)))
if 1:
    mgo.calc = GPAW(mode=PW(500),
                    parallel=dict(band=1),
                    idiotproof=False,
                    communicator=comm,
                    setups={'Mg': '2'},
                    convergence={'eigenstates': 5.e-9},
                    kpts=monkhorst_pack((2, 2, 2)) + 0.25)
    mgo.get_potential_energy()
    if rank < 3:
       mgo.calc.write('mgo', 'all')
    else:
       mgo.calc.write('dummy_%d' % rank, 'all')
    world.barrier()

for name in ['PBE0', 'HSE03', 'HSE06']:
    calc = GPAW('mgo', setups={'Mg': '2'},
               txt=None, communicator=serial_comm)
    hyb_calc = HybridXC(name, alpha=5.0, bandstructure=True, world=comm)
    de_skn = vxc(calc, hyb_calc) - vxc(calc, 'LDA')
        
    if name == 'PBE0':
        de_skn_test = np.array([-1.3700, -1.3643, -1.3777, -24.184590])
    if name == 'HSE03':
        de_skn_test = np.array([-2.4565, -2.4326, -2.4583, -24.405001])
    if name == 'HSE06':
        de_skn_test = np.array([-2.0311, -2.0151, -2.0367, -24.324485])

    if rank == 0:
        print(de_skn[0, 0, 1:4], abs(de_skn[0, 0, 1:4] - de_skn_test[0]).max())
        print(de_skn[0, 1, 2:4], abs(de_skn[0, 1, 2:4] - de_skn_test[1]).max())
        print(de_skn[0, 2, 2:4], abs(de_skn[0, 2, 2:4] - de_skn_test[2]).max())
        print(hyb_calc.exx, abs(hyb_calc.exx - de_skn_test[3]))
        assert abs(de_skn[0, 0, 1:4] - de_skn_test[0]).max() < 0.02
        assert abs(de_skn[0, 1, 2:4] - de_skn_test[1]).max() < 0.008
        assert abs(de_skn[0, 2, 2:4] - de_skn_test[2]).max() < 0.004
        assert abs(hyb_calc.exx - de_skn_test[3]) < 2e-4
