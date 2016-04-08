from __future__ import print_function
import numpy as np
from ase.lattice import bulk
from ase.dft.kpoints import monkhorst_pack
from ase.units import Bohr

from gpaw import GPAW
from gpaw.response.chi import CHI
from gpaw.response.chi0 import Chi0
from gpaw.mpi import serial_comm


omega = np.array([0, 1.0, 2.0])
for k in [2, 3]:
    q_c = [0, 0, 1.0 / k]
    for gamma in [False, True]:
        if k == 3 and gamma:
            continue
        kpts = monkhorst_pack((k, k, k))
        if gamma:
            kpts += 0.5 / k
        for center in [False, True]:
            a = bulk('Si', 'diamond')
            if center:
                a.center()
            for sym in [False, True]:
                name = 'si.k%d.g%d.c%d.s%d' % (k, gamma, center, bool(sym))
                print(name)
                if 1:
                    calc = a.calc = GPAW(
                        kpts=kpts,
                        eigensolver='rmm-diis',
                        symmetry={'point_group': sym},
                        mode='pw',
                        width=0.001,
                        txt=name + '.txt')
                    e = a.get_potential_energy()
                    #calc.diagonalize_full_hamiltonian(nbands=100)
                    calc.write(name, 'all')
                    
                calc = GPAW(name, txt=None, communicator=serial_comm)

                chiold = CHI(calc, w=omega, q=q_c, ecut=100,
                             hilbert_trans=False, xc='RPA',
                             G_plus_q=True, txt=name + '.logold')
                chiold.initialize()
                chiold.calculate()
                chi0old_wGG = chiold.chi0_wGG

                chi = Chi0(calc, omega, hilbert=False,
                           ecut=100, txt=name + '.log')
                pd, chi0_wGG, _, _ = chi.calculate(q_c)

                assert abs(chi0_wGG - chi0old_wGG).max() < 1e-15
                
                if not sym and not center:
                    chi00_w = chi0_wGG[:, 0, 0]
                elif -1 not in calc.wfs.kd.bz2bz_ks:
                    assert abs(chi0_wGG[:, 0, 0] - chi00_w).max() < 3e-5
                    #print abs(chi0_wGG[:, 0, 0] - chi00_w).max()
                    
                if not sym:
                    chi00_wGG = chi0_wGG
                elif -1 not in calc.wfs.kd.bz2bz_ks:
                    assert abs(chi0_wGG - chi00_wGG).max() < 2e-5
                    #print abs(chi0_wGG - chi00_wGG).max()

                q0_c = [0, 1e-7, 1e-7]
                q0_v = np.dot(q0_c, a.get_reciprocal_cell() * 2 * np.pi) * Bohr
                q0 = (q0_v**2).sum()**0.5
                
                chiold = CHI(calc, w=omega, q=q0_c, ecut=100,
                             hilbert_trans=False, xc='RPA', optical_limit=True,
                             G_plus_q=True, txt=name + '.logold0')
                chiold.initialize()
                chiold.calculate()
                chi0old_wGG = chiold.chi0_wGG
                chi0old_wGG[:, 0] /= q0
                chi0old_wGG[:, :, 0] /= q0
                
                pd, chi0_wGG, _, _ = chi.calculate([0, 0, 0])
                
                assert abs(chi0_wGG - chi0old_wGG).max() < 0.003
                assert abs(chi0_wGG - chi0old_wGG)[:, 1:, 1:].max() < 1e-9
                
                if not sym and not center:
                    chi000_w = chi0_wGG[:, 0, 0]
                elif -1 not in calc.wfs.kd.bz2bz_ks:
                    assert abs(chi0_wGG[:, 0, 0] - chi000_w).max() < 0.001
                    
                if not sym:
                    chi000_wGG = chi0_wGG
                elif -1 not in calc.wfs.kd.bz2bz_ks:
                    assert abs(chi0_wGG - chi000_wGG).max() < 0.001
