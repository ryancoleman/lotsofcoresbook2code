from ase import Atoms
from gpaw import GPAW
from gpaw.xc.vdw import VDWFunctional

vdw = VDWFunctional('vdW-DF', verbose=1)
L = 2.5
a = Atoms('H', cell=(L, L, L), pbc=True, calculator=GPAW(nbands=1))
e = a.get_potential_energy()
e2 = a.calc.get_xc_difference(vdw)
print(e, e2)
assert (vdw.shape == (24, 24, 24)).all()
