from __future__ import print_function
import ase.db
from ase.units import kcal, mol
from ase.data.g2_1_ref import ex_atomization, atomization


c = ase.db.connect('results.db')

# Energy of atoms:
atoms = {}
for d in c.select(natoms=1):
    atoms[d.numbers[0]] = [d.energy, d.exx]
    
maepbe = 0.0
maeexx = 0.0
print('                 PBE                   EXX')
print(('-' * 48))
for d in c.select('natoms>1'):
    epberef = atomization[d.name][2] * kcal / mol
    eexxref = ex_atomization[d.name][0] * kcal / mol

    epbe = sum(atoms[atom][0] for atom in d.numbers) - d.energy
    eexx = sum(atoms[atom][1] for atom in d.numbers) - d.exx

    maepbe += abs(epbe - epberef) / len(ex_atomization)
    maeexx += abs(eexx - eexxref) / len(ex_atomization)

    print(('%-4s %10.3f %10.3f %10.3f %10.3f' %
          (d.name, epbe, epbe - epberef, eexx, eexx - eexxref)))

print(('-' * 48))
print(('MAE  %10.3f %10.3f %10.3f %10.3f' % (0, maepbe, 0, maeexx)))

assert maepbe < 0.025
assert maeexx < 0.05
