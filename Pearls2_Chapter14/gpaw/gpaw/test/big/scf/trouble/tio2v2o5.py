from ase.io import read
from gpaw import GPAW
atoms = read('tio2v2o5.out')
atoms.calc = GPAW(h=0.20,
                  kpts=(2,2,1),
                  xc='RPBE')
ncpus = 8
