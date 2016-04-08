import numpy as np
import cPickle as pickle
from ase.structure import molecule
from gpaw.lcao.tools import makeU, makeV
from gpaw import GPAW, FermiDirac, restart
from gpaw.lcao.pwf2 import LCAOwrap
from gpaw.lcao.tools import remove_pbc, get_bfi2, get_bf_centers
from gpaw.mpi import world, rank, MASTER, serial_comm
from gpaw.test import equal

scat = range(2)
atoms = molecule('H2')
atoms.set_cell([6.4, 6.4, 6.4])
atoms.center()
calc = GPAW(mode='lcao', occupations=FermiDirac(0.1))
atoms.set_calculator(calc)
atoms.get_potential_energy()
calc.write('lcao_pair.gpw')

if rank == MASTER:
    atoms, calc = restart('lcao_pair.gpw',
                          txt=None, communicator=serial_comm)
    lcao = LCAOwrap(calc)
    fermi = calc.get_fermi_level()
    H = lcao.get_hamiltonian()
    S = lcao.get_overlap()
    pickle.dump((H, S), open('lcao_pair_hs.pckl', 'wb'), 2)
    symbols = atoms.get_chemical_symbols()
    #indices = get_bfi2(symbols, basis, scat)
    indices = range(2)
    lcao.get_xc(indices=indices).dump('lcao_pair_xc.pckl')
    lcao.get_Fcore(indices=indices).dump('lcao_pair_Fcore.pckl')
    w_wG = lcao.get_orbitals(indices=indices)
    P_awi = lcao.get_projections(indices=indices)
    pickle.dump((w_wG, P_awi), open('lcao_pair_w_wG__P_awi.pckl', 'wb'), 2)

world.barrier()
makeU('lcao_pair.gpw',
      'lcao_pair_w_wG__P_awi.pckl',
      'lcao_pair_eps_q__U_pq.pckl',
      1.0e-5,
      dppname='lcao_pair_D_pp.pckl')

world.barrier()
makeV('lcao_pair.gpw',
      'lcao_pair_w_wG__P_awi.pckl',
      'lcao_pair_eps_q__U_pq.pckl',
      'lcao_pair_V_qq.pckl',
      'lcao_pair_V_qq.log',
      False)

world.barrier()
V_qq = np.load('lcao_pair_V_qq.pckl')
eps_q, U_pq = np.load('lcao_pair_eps_q__U_pq.pckl')
assert U_pq.flags.contiguous
Usq_pq = U_pq * np.sqrt(eps_q)
V_pp = np.dot(np.dot(Usq_pq, V_qq), Usq_pq.T.conj())
V_pp_ref = np.array(
[[ 15.34450177,  11.12669608,  11.12669608,  12.82934563],
 [ 11.12669608,   8.82280293,   8.82280293,  11.12669608],
 [ 11.12669608,   8.82280293,   8.82280293,  11.12669608],
 [ 12.82934563,  11.12669608,  11.12669608,  15.34450178]])
equal(abs(V_pp_ref-V_pp).max(), 0.0, 1.0e-5)
