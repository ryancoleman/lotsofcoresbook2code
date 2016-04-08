from gpaw import GPAW
from ase.lattice.surface import fcc100


def energy(N, k, a=4.05):
    fcc = fcc100('Al', (1, 1, N), a=a, vacuum=7.5)
    fcc.center(axis=2)
    calc = GPAW(nbands=N * 3,
                kpts=(k, k, 1),
                h=0.25,
                txt='slab-%d.txt' % N)
    fcc.set_calculator(calc)
    e = fcc.get_potential_energy()
    calc.write('slab-%d.gpw' % N)
    return e

e = energy(N, k)
