from ase import Atoms
from gpaw import GPAW, FermiDirac

L = 2.5
a = Atoms('H', cell=(L, L, L), pbc=True)
calc = GPAW(xc='vdW-DF',
            occupations=FermiDirac(width=0.001),
            txt='H.vdW-DF.txt')
a.set_calculator(calc)
e1 = a.get_potential_energy()

calc.set(txt='H.vdW-DF.spinpol.txt',
         spinpol=True,
         occupations=FermiDirac(width=0.001, fixmagmom=True))
e2 = a.get_potential_energy()

assert abs(calc.get_eigenvalues(spin=0)[0] -
           calc.get_eigenvalues(spin=1)[0]) < 1e-10

assert abs(e1 - e2) < 5e-8, abs(e1 - e2)
