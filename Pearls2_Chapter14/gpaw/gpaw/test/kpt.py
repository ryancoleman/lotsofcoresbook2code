import numpy as np
from ase.dft.kpoints import monkhorst_pack
from gpaw.kpt_descriptor import KPointDescriptor, to1bz
k = 70
k_kc = monkhorst_pack((k, k, 1))
kd = KPointDescriptor(k_kc + (0.5 / k, 0.5 / k, 0))
assert (kd.N_c == (k, k, 1)).all()
assert abs(kd.offset_c - (0.5 / k, 0.5 / k, 0)).sum() < 1e-9

bzk_kc = np.array([[0.5, 0.5, 0],
                   [0.50000000001, 0.5, 0],
                   [0.49999999999, 0.5, 0],
                   [0.55, -0.275, 0]])
cell_cv = np.array([[1, 0, 0],
                    [-0.5, 3**0.5 / 2, 0],
                    [0, 0, 5]])
bz1k_kc = to1bz(bzk_kc, cell_cv)
error_kc = bz1k_kc - np.array([[0.5, -0.5, 0],
                               [0.50000000001, -0.5, 0],
                               [0.49999999999, -0.5, 0],
                               [0.55, -0.275, 0]])
assert abs(error_kc).max() == 0.0
assert KPointDescriptor(np.zeros((1, 3)) + 1e-14).gamma
