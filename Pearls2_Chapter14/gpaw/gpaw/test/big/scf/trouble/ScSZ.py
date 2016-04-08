from ase.io import read
from gpaw import GPAW, ConvergenceError
from gpaw.mixer import Mixer
from gpaw.utilities import compiled_with_sl

# the system loosing magnetic moment

atoms = read('ScSZ.xyz')
atoms.set_cell([[  7.307241,   0.,         0.,      ],
               [  0.,        12.656514,   0.,      ],
               [  0.,         0.,        19.,      ]],
              scale_atoms=False)
atoms.center(axis=2)
magmoms = [0.0 for n in range(len(atoms))]
for n, a in enumerate(atoms):
    if a.symbol == 'Ni':
        magmoms[n] = 0.6
atoms.set_initial_magnetic_moments(magmoms)

atoms.pbc = (True, True, False)

atoms.calc = GPAW(h=0.20,
                  kpts=(2,1,1),
                  xc='PBE',
                  width=0.1,
                  maxiter=100,
                  txt='ScSZ.txt',
                  )

ncpus = 8
