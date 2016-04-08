from ase import Atoms
from gpaw import GPAW

# Benzene ring:
benzene = Atoms(symbols='C6H6',
                positions=[
    ( 0.000000,  1.395248, 0.000000),
    ( 1.208320,  0.697624, 0.000000),
    ( 1.208320, -0.697624, 0.000000),
    ( 0.000000, -1.395248, 0.000000),
    (-1.208320, -0.697624, 0.000000),
    (-1.208320,  0.697624, 0.000000),
    ( 0.000000,  2.482360, 0.000000),
    ( 2.149787,  1.241180, 0.000000),
    ( 2.149787, -1.241180, 0.000000),
    ( 0.000000, -2.482360, 0.000000),
    (-2.149787, -1.241180, 0.000000),
    (-2.149787,  1.241180, 0.000000)])

benzene.center(vacuum=2.5)

calc = GPAW(nbands=15,
            h=0.2,
            xc='PBE',
            txt='benzene.txt')

benzene.set_calculator(calc)
benzene.get_potential_energy()
calc.write('benzene.gpw', mode='all')
