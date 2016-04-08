from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
def groundstate(a, k):
    si = bulk('Si', 'diamond', a)
    si.calc = GPAW(kpts=(k, k, k),
                   xc='PBE',
                   gpts=(20, 20, 20),
                   occupations=FermiDirac(0.01),
                   parallel={'domain': 1, 'band': 1},
                   txt='Si-PBE-%.3f-%d.txt' % (a, k))
    si.get_potential_energy()
    return si
