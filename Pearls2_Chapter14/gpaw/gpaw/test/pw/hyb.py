from __future__ import print_function
from ase import Atoms
from gpaw import GPAW, PW
from gpaw.mpi import rank, size, serial_comm
from gpaw.xc.hybridg import HybridXC

a = 2.0
li = Atoms('Li', cell=(a, a, a), pbc=1)
for spinpol in [False, True]:
    for symm in [{}, 'off', {'time_reversal': False, 'point_group': False}]:
        if size == 8 and not spinpol and symm == {}:
            continue
        for qparallel in [False, True]:
            if rank == 0:
                print((spinpol, symm, qparallel))
            li.calc = GPAW(mode=PW(300),
                           kpts=(2, 3, 4),
                           spinpol=spinpol,
                           symmetry=symm,
                           parallel={'band': 1},
                           txt=None,
                           idiotproof=False)
            e = li.get_potential_energy()
            if qparallel:
                li.calc.write('li', mode='all')
                calc = GPAW('li', txt=None, communicator=serial_comm)
            else:
                calc = li.calc
            exx = HybridXC('EXX',
                           logfilename=None,
                           method='acdf')
            de = calc.get_xc_difference(exx)
            exx = HybridXC('EXX',
                           logfilename=None,
                           method='acdf',
                           bandstructure=True, bands=[0, 1])
            de2 = calc.get_xc_difference(exx)
            kd = calc.wfs.kd
            print(e, -0.56024, abs(e - -0.56024))
            print(de, -0.4520, abs(de - -0.4520))
            print(de, de2, abs(de - de2))
            assert abs(e - -0.56024) < 1e-5, abs(e)
            assert abs(de - -0.4520) < 3e-4, abs(de)
            assert abs(de - de2) < 1e-12
            for k in range(kd.nibzkpts):
                if abs(kd.ibzk_kc[k] - [0.25, 1 / 3.0, 3 / 8.0]).max() < 1e-7:
                    assert abs(exx.exx_skn[:, k, 0] - -0.18246).max() < 1e-5
