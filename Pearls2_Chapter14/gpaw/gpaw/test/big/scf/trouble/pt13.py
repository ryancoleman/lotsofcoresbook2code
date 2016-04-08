from ase.lattice.cubic import FaceCenteredCubic
from gpaw import GPAW, PW
element = 'Pt'
atoms = FaceCenteredCubic(symbol=element,
                          size=(2, 2, 2),
                          directions=[[1, 1, 0],
                                      [-1, 1, 0],
                                      [0, 0, 1]])
del atoms[4]
del atoms[3]
del atoms[2]

h = 0.16
kpts=(8, 8, 4)
ecut = 800
xc1 = 'PBE'
atoms.calc = GPAW(mode=PW(ecut=ecut),
                  kpts=kpts,
                  xc=xc1)
ncpus = 8
