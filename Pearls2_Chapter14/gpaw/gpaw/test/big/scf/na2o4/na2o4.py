import numpy as np

from ase.io import read
from ase.dft.kpoints import monkhorst_pack

from gpaw import GPAW, ConvergenceError
from gpaw.wavefunctions.pw import PW
from gpaw.mixer import Mixer

# the system does not converge with 3 diis steps, but converges with 2, 4, ...

bulk = read('na2o4.xyz')
bulk.set_cell([[ 4.31110615,  0.,          0.        ],
               [ 0.,          5.60646198,  0.        ],
               [ 0.,          0.,          3.48126881]],
              scale_atoms=False)
magmoms = [0.0 for n in range(len(bulk))]
for n, a in enumerate(bulk):
    if a.symbol == 'O':
        magmoms[n] = 0.5
bulk.set_initial_magnetic_moments(magmoms)

bulk.pbc = (True, True, True)

nk = 4
kpts = monkhorst_pack((nk,nk,nk))
kshift = 1./(2*nk)
kpts += np.array([kshift, kshift, kshift])

calc = GPAW(
    mode=PW(),
    xc='PBE',
    kpts=kpts,
    parallel={'band': 1},
    eigensolver='cg',
    txt='na2o4.txt',
    )

bulk.set_calculator(calc)

try:
    bulk.get_potential_energy()
except ConvergenceError:
    pass

assert not calc.scf.converged

del bulk.calc

calc1 = GPAW(
    mode=PW(),
    xc='PBE',
    kpts=kpts,
    parallel={'band': 1},
    eigensolver='cg',
    txt='na2o4.txt',
    )

calc1.set(mixer=Mixer(0.10,2))
calc1.set(txt='na2o4_m1.txt')
bulk.set_calculator(calc1)
bulk.get_potential_energy()
