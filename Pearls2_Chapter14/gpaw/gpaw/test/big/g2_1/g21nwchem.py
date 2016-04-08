import ase.db
from ase.structure import molecule
from ase.optimize.bfgs import BFGS
from ase.calculators.nwchem import NWChem
from ase.data.g2_1 import molecule_names, atom_names


c = ase.db.connect('g2-1.db')

for name in molecule_names + atom_names:
    id = c.reserve(name=name, calculator='nwchem')
    if id is None:
        continue
    atoms = molecule(name)
    atoms.calc = NWChem(command='mpiexec -np 4 nwchem PREFIX.nw > PREFIX.out',
                        geometry='noautosym nocenter noautoz',
                        task='gradient',
                        xc='PBE',
                        grid='nodisk',
                        tolerances='tight',
                        basis='def2-qzvppd',
                        basispar='spherical',
                        direct='noio',
                        label=name)
    atoms.get_forces()
    c.write(atoms, name=name, relaxed=False)
    if len(atoms) > 1:
        opt = BFGS(atoms, logfile=name + '.nwchem.log')
        opt.run(0.01)
        c.write(atoms, name=name, relaxed=True)
    del c[id]
