from ase.io import read
from gpaw import GPAW, FermiDirac
ncpus = 8
atoms = read('biimtf.xyz')
atoms.center(vacuum=5)
atoms.calc = GPAW(h=0.16,
                  charge=+1,
                  occupations=FermiDirac(0.05),
                  xc='RPBE')
