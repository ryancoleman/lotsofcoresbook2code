from math import log
from ase import Atoms
from ase.units import Bohr
from gpaw import GPAW, FermiDirac
from gpaw.test import equal

a = 4.0
h = 0.2
hydrogen = Atoms('H',
                 [(a / 2, a / 2, a / 2)],
                 cell=(a, a, a))

hydrogen.calc = GPAW(h=h, nbands=1, convergence={'energy': 1e-7})
e1 = hydrogen.get_potential_energy()
equal(e1, 0.526939, 0.001)

dens = hydrogen.calc.density
c = dens.gd.find_center(dens.nt_sG[0]) * Bohr
equal(abs(c - a / 2).max(), 0, 1e-13)

kT = 0.001
hydrogen.calc.set(occupations=FermiDirac(width=kT))
e2 = hydrogen.get_potential_energy()
equal(e1, e2 + log(2) * kT, 3.0e-7)
