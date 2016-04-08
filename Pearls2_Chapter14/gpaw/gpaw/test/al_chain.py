from ase import Atoms
from gpaw import GPAW
from gpaw.mpi import size
d = 4.0 / 2**0.5
ndomains = size // 8 + 1
calc = GPAW(h=d / 16, kpts=(17, 1, 1), parallel={'domain': ndomains,
                                                 'band': 1})
chain = Atoms('Al', cell=(d, 5, 5), pbc=True, calculator=calc)
e = chain.get_potential_energy()
assert abs(e - -1.81816) < 0.00005
assert calc.wfs.kd.comm.size * ndomains == size
calc.write('al_chain', 'all')
