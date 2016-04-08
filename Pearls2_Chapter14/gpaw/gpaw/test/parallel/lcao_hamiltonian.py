import cPickle as pickle
import numpy as np
from ase import Atoms
from gpaw import GPAW, restart, setup_paths
from gpaw.lcao.tools import get_lcao_hamiltonian
from gpaw.mpi import world
from gpaw.atom.basis import BasisMaker
from gpaw.test import equal

if world.rank == 0:
    basis = BasisMaker('Li', 'szp').generate(1, 1)
    basis.write_xml()
world.barrier()
if setup_paths[0] != '.':
    setup_paths.insert(0, '.')

if 1:
    a = 2.7
    bulk = Atoms('Li', pbc=True, cell=[a, a, a])
    calc = GPAW(gpts=(8, 8, 8), kpts=(4, 4, 4), mode='lcao', basis='szp')
    bulk.set_calculator(calc)
    e = bulk.get_potential_energy()
    niter = calc.get_number_of_iterations()
    calc.write('temp.gpw')

atoms, calc = restart('temp.gpw')
H_skMM, S_kMM = get_lcao_hamiltonian(calc)
eigs = calc.get_eigenvalues(kpt=2)

if world.rank == 0:
    eigs2 = np.linalg.eigvals(np.linalg.solve(S_kMM[2], H_skMM[0, 2])).real
    eigs2.sort()
    assert abs(sum(eigs - eigs2)) < 1e-8

    energy_tolerance = 0.00003
    niter_tolerance = 0
    equal(e, -1.82847, energy_tolerance)
    equal(niter, 5, niter_tolerance)
