import ase.db
from ase.structure import molecule
from ase.optimize.bfgs import BFGS
from ase.data.g2_1 import molecule_names, atom_names

from gpaw import GPAW, PW, Mixer


c = ase.db.connect('g2-1.db')

for name in molecule_names + atom_names:
    id = c.reserve(name=name, calculator='gpaw')
    if id is None:
        continue
    atoms = molecule(name)
    atoms.cell = [12, 13.01, 14.02]
    atoms.center()
    atoms.calc = GPAW(mode=PW(800),
                      xc='PBE',
                      fixmom=True,
                      txt=name + '.txt')
    if name in ['Na2', 'NaCl', 'NO', 'ClO', 'Cl', 'Si']:
        atoms.calc.set(eigensolver='dav',
                       mixer=Mixer(0.05, 2))
    atoms.get_forces()
    c.write(atoms, name=name, relaxed=False)
    if len(atoms) > 1:
        opt = BFGS(atoms, logfile=name + '.gpaw.log')
        opt.run(0.01)
        c.write(atoms, name=name, relaxed=True)
    del c[id]
