import numpy as np
from ase.lattice import bulk
from ase.units import Hartree, Bohr
from gpaw import GPAW, FermiDirac
from ase.dft.kpoints import monkhorst_pack


for kpt in (3, 4):
    kpts = (kpt, kpt, kpt)
    bzk_kc = monkhorst_pack(kpts)
    shift_c = []
    for Nk in kpts:
        if Nk % 2 == 0:
            shift_c.append(0.5 / Nk)
        else:
            shift_c.append(0.)
    bzk_kc += shift_c

    atoms = bulk('Si', 'diamond', a=5.431)
    calc = GPAW(h=0.2,
                kpts=bzk_kc)
    
    atoms.set_calculator(calc)
    atoms.get_potential_energy()

    kd = calc.wfs.kd
    bzq_qc = kd.get_bz_q_points()
    ibzq_qc = kd.get_ibz_q_points(bzq_qc, calc.wfs.kd.symmetry.op_scc)[0]

    assert np.abs(bzq_qc - kd.bzk_kc).sum() < 1e-8
    assert np.abs(ibzq_qc - kd.ibzk_kc).sum() < 1e-8
