from __future__ import print_function
from ase import Atoms
from gpaw import GPAW
from gpaw.utilities.sic import NSCFSIC

atoms = ['He','Be'] #,'Ne'] # Ne deviates already 2.5 eV
EE = []
EREF = [-79.4,-399.8,-3517.6]

for a in atoms:
    s = Atoms(a)
    s.center(vacuum=4.0)
    calc = GPAW(h=0.15, txt=a + '.txt')
    s.set_calculator(calc)
    E = s.get_potential_energy()
    EE.append(NSCFSIC(calc).calculate())

print("Difference to table VI of Phys. Rev. B 23, 5048 in eV")
#http://prola.aps.org/abstract/PRB/v23/i10/p5048_1
print("%10s%10s%10s%10s" % ("atom", "ref.", "gpaw", "diff"))
for a, er, e in zip(atoms, EREF, EE):
    print("%10s%10.2f%10.2f%10.2f" % (a, er, e, er-e))
    assert abs(er-e)<0.1 
    # Arbitary 0.1 eV tolerance for non-self consistent SIC
    # Note that Ne already deviates 2.5 eV
