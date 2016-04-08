from ase import Atoms
from gpaw import GPAW, PW
h = Atoms('H', cell=(5, 5, 5))
h.center()
h.calc = GPAW(setups='ae', txt='H.ae.txt')
for ecut in range(200, 1001, 100):
    h.calc.set(mode=PW(ecut))
    e = h.get_potential_energy()
