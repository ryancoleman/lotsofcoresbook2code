"""Test band-gaps for Si."""

from ase.lattice import bulk
import numpy as np

from gpaw.test import equal
from gpaw import GPAW, PW, FermiDirac
from gpaw.response.g0w0 import G0W0


def run(atoms, symm, name):
    atoms.calc = GPAW(mode=PW(250),
                      eigensolver='rmm-diis',
                      occupations=FermiDirac(0.01),
                      symmetry=symm,
                      kpts={'size': (2, 2, 2), 'gamma': True},
                      txt=name + '.txt')
    e = atoms.get_potential_energy()
    scalapack = atoms.calc.wfs.bd.comm.size
    atoms.calc.diagonalize_full_hamiltonian(nbands=8, scalapack=scalapack)
    atoms.calc.write(name, mode='all')
    gw = G0W0(name, 'gw-' + name,
              nbands=8,
              kpts=[(0, 0, 0), (0.5, 0.5, 0)],  # Gamma, X
              ecut=40,
              domega0=0.1,
              eta=0.2,
              bands=(3, 7),  # homo, lumo, lumo+1, lumo+2
              )
    results = gw.calculate()
    return e, results

    
a = 5.43
si1 = bulk('Si', 'diamond', a=a)
si2 = si1.copy()
si2.positions -= a / 8

i = 0
results = []
for si in [si1, si2]:
    for symm in [{}, 'off', {'time_reversal': False}, {'point_group': False}]:
        e, r = run(si, symm, str(i))
        G, X = r['eps'][0]
        results.append([e, G[0], G[1] - G[0], X[1] - G[0], X[2] - X[1]])
        G, X = r['qp'][0]
        results[-1].extend([G[0], G[1] - G[0], X[1] - G[0], X[2] - X[1]])
        i += 1

equal(abs(np.array(results[0]) -
          [-9.25,
           5.44, 2.39, 0.40, 0,
           6.26, 3.57, 1.32, 0]).max(), 0, 0.01)
equal(np.ptp(results, 0).max(), 0, 0.005)
