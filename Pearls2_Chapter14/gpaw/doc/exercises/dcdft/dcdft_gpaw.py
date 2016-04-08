from time import time
import numpy as np
import ase.db
from ase.test.tasks.dcdft import DeltaCodesDFTCollection as Collection
from gpaw import GPAW, PW, FermiDirac

c = ase.db.connect('dcdft.db')

ecut = 340
kptdensity = 3.5
width = 0.10

collection = Collection()

for name in ['K', 'Ca', 'Ti']:
    atoms = collection[name]
    cell = atoms.get_cell()
    # Loop over volumes:
    for n, x in enumerate(np.linspace(0.98, 1.02, 5)):
        id = c.reserve(name=name, x=x)
        if id is None:
            # This calculation has been or is being done:
            continue

        atoms.set_cell(cell * x, scale_atoms=True)
        atoms.calc = GPAW(txt='%s-%d.txt' % (name, n),
                          mode=PW(ecut),
                          xc='PBE',
                          kpts={'density': kptdensity},
                          occupations=FermiDirac(width))

        t1 = time()
        atoms.get_potential_energy()
        t2 = time()

        # Write to database:
        c.write(atoms, name=name, x=x, time=t2 - t1,
                ecut=ecut, kptdensity=kptdensity, width=width)

        del c[id]
