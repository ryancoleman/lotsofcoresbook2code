from gpaw import GPAW, PW, FermiDirac
from ase.lattice import bulk
si = bulk('Si')
si.calc = GPAW(mode=PW(200),
               kpts={'size': (2, 2, 2), 'gamma': True},
               occupations=FermiDirac(0.01))
si.get_potential_energy()
si.calc.diagonalize_full_hamiltonian()
si.calc.write('Si_gs', 'all')
dct = {}
execfile('Si_g0w0_ppa.py', dct)
assert abs(dct['ks_gap'] - 0.404) < 0.01
assert abs(dct['qp_gap'] - 1.138) < 0.01
