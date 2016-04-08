from math import sqrt
import numpy as np

from gpaw.symmetry import Symmetry
from ase.dft.kpoints import monkhorst_pack

# Primitive diamond lattice, with Si lattice parameter
a = 5.475
cell_cv = .5 * a * np.array([(1, 1, 0), (1, 0, 1), (0, 1, 1)])
spos_ac = np.array([(.00, .00, .00),
                    (.25, .25, .25)])
id_a = [1, 1] # Two identical atoms
pbc_c = np.ones(3, bool)
bzk_kc = monkhorst_pack((4, 4, 4))

# Do check
symm = Symmetry(id_a, cell_cv, pbc_c, symmorphic=False)
symm.analyze(spos_ac)
ibzk_kc, w_k = symm.reduce(bzk_kc)[:2]
assert len(symm.op_scc) == 48
assert len(w_k) == 10
a = 3 / 32.; b = 1 / 32.; c = 6 / 32.
assert np.all(w_k == [a, b, a, c, c, a, a, a, a, b])
assert not symm.op_scc.sum(0).any()

# Rotate unit cell and check again:
cell_cv = a / sqrt(2) * np.array([(1, 0, 0),
                                  (0.5, sqrt(3) / 2, 0),
                                  (0.5, sqrt(3) / 6, sqrt(2.0 / 3))])
symm = Symmetry(id_a, cell_cv, pbc_c, symmorphic=False)
symm.analyze(spos_ac)
ibzkb_kc, wb_k = symm.reduce(bzk_kc)[:2]
assert len(symm.op_scc) == 48
assert abs(w_k - wb_k).sum() < 1e-14
assert abs(ibzk_kc - ibzkb_kc).sum() < 1e-14
assert not symm.op_scc.sum(0).any()

bzk_kc = monkhorst_pack((3, 3, 3))
symm = Symmetry(id_a, cell_cv, pbc_c)
symm.analyze(spos_ac)
ibzk_kc, w_k = symm.reduce(bzk_kc)[:2]
assert len(symm.op_scc) == 24
assert len(w_k) == 4
assert abs(w_k * 27 - (1, 12, 6, 8)).sum() < 1e-14
assert not symm.op_scc.sum(0).any()

# Rocksalt Ni2O2
a = 7.92; x = 2. * np.sqrt(1./3.); y = np.sqrt(1./8.); z1 = np.sqrt(1./24.); z2 = np.sqrt(1./6.)
cell_cv = a * np.array([(x, y, -z1), (x, -y, -z1), (x, 0., z2)])
spos_ac = np.array([[0., 0. ,0.], [1./2., 1./2., 1./2.], [1./4., 1./4., 1./4.], [3./4., 3./4., 3./4.]])
id_a = [1, 2, 3, 3]
pbc_c = np.array([1, 1, 1], bool)
bzk_kc = monkhorst_pack((2, 2, 2))

# Do check
symm = Symmetry(id_a, cell_cv, pbc_c, symmorphic=False)
symm.analyze(spos_ac)
ibzk_kc, w_k = symm.reduce(bzk_kc)[:2]
assert len(symm.op_scc) == 12
assert len(w_k) == 2
assert np.all(w_k == [3/4., 1/4.])

#AlF3
a = 2.465250000; b = 1.423312751; c = 4.148733333
a1 = np.array([     a,-1.0*b, c]); a2 = np.array([ 0.0  , 2.0*b, c]); a3 = np.array([-1.0*a,-1.0*b, c])
cell_cv = np.array([a1,a2,a3])
id_a = [1, 1, 2, 2, 2, 2, 2, 2]
spos_ac = np.array([(0.000000000,   0.000000000,   0.000000000), (0.500000000,   0.500000000,   0.500000000),
                    (0.647758346,   0.852241654,   0.250000000), (0.852241654,   0.250000000,   0.647758346),
                    (0.250000000,   0.647758346,   0.852241654), (0.352241654,   0.147758346,   0.750000000),
                    (0.147758346,   0.750000000,   0.352241654), (0.750000000,   0.352241654,   0.147758346)])
pbc_c = np.array([1, 1, 1], bool)
bzk_kc = monkhorst_pack((3, 3, 3))

# Do check
symm = Symmetry(id_a, cell_cv, pbc_c, symmorphic=False)
symm.analyze(spos_ac)
ibzk_kc, w_k = symm.reduce(bzk_kc)[:2]
assert len(symm.op_scc) == 12
assert len(w_k) == 6
assert np.all(w_k == [1/27., 2/9., 2/9., 2/9., 2/9., 2/27.])
assert not symm.op_scc.sum(0).any()
