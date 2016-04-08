from __future__ import print_function
import numpy as np
from ase.structure import molecule

from gpaw import GPAW
from gpaw.eigensolvers.rmm_diis_old import RMM_DIIS
from gpaw.lcao.projected_wannier import get_lcao_projections_HSP

atoms = molecule('C2H2')
atoms.center(vacuum=3.0)
calc = GPAW(gpts=(32, 32, 48), eigensolver=RMM_DIIS())
atoms.set_calculator(calc)
atoms.get_potential_energy()

V_qnM, H_qMM, S_qMM, P_aqMi = get_lcao_projections_HSP(
    calc, bfs=None, spin=0, projectionsonly=False)


# Test H and S
eig = np.linalg.eigvals(np.linalg.solve(S_qMM[0], H_qMM[0])).real
eig.sort()
eig_ref = np.array([-17.87911292, -13.24864985, -11.43106707, -7.12558127,
                    -7.12558127, 0.59294531, 0.59294531, 3.92526888,
                    7.45117399, 26.73466374])
print(eig)
assert np.allclose(eig, eig_ref)

# Test V
Vref_nM = np.array(
    [[-0.845  , -0.     ,  0.39836,  0.     , -0.845  ,  0.     ,
      -0.39836, -0.     , -0.40865, -0.40865],
     [-0.43521,  0.     , -0.61583, -0.     ,  0.43521,  0.     ,
      -0.61583, -0.     ,  0.60125, -0.60125],
     [ 0.12388,  0.     ,  0.62178,  0.     ,  0.12388,  0.     ,
      -0.62178,  0.     ,  0.50672,  0.50672],
     [-0.     , -0.56118,  0.     ,  0.62074, -0.     , -0.56118,
      -0.     ,  0.62074, -0.     , -0.     ],
     [-0.     ,  0.62074, -0.     ,  0.56118,  0.     ,  0.62074,
       0.     ,  0.56118, -0.     , -0.     ],
     [ 0.     ,  0.5003 ,  0.     ,  0.10456,  0.     , -0.5003 ,
      -0.     , -0.10456, -0.     , -0.     ],
     [-0.     ,  0.10456,  0.     , -0.5003 , -0.     , -0.10456,
      -0.     ,  0.5003 ,  0.     ,  0.     ],
     [ 0.12449, -0.     , -0.02292,  0.     , -0.12449,  0.     ,
      -0.02292, -0.     ,  0.29925, -0.29925],
     [-0.05659,  0.     , -0.02315,  0.     , -0.05659, -0.     ,
       0.02315, -0.     ,  0.27325,  0.27325],
     [-0.15622,  0.     ,  0.14149, -0.     ,  0.15622,  0.     ,
       0.14149, -0.     ,  0.04507, -0.04507]])
## np.set_printoptions(precision=5, suppress=1)
## from gpaw.mpi import rank
## if rank == 0:
##     print V_qnM[0]
print(abs(V_qnM[0]) - abs(Vref_nM)) # assert <--- this to zero
