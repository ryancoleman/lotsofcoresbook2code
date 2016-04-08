from ase import Atoms
from gpaw import GPAW, PW, FermiDirac

# Setup up bulk NiO in an antiferromagnetic configuration:
a = 4.19  # lattice constants
b = a / 2**0.5
m = 2.0
atoms = Atoms('Ni2O2',
              pbc=True,
              cell=(b, b, a),
              positions=[(0, 0, 0),
                         (b / 2, b / 2, a / 2),
                         (0, 0, a / 2),
                         (b / 2, b / 2, 0)],
              magmoms=[m, -m, 0, 0])

k = 2  # number of k-points
atoms.calc = GPAW(mode=PW(400),
                  occupations=FermiDirac(width=0.05),
                  setups={'Ni': ':d,6.0'},  # U=6 eV for Ni d orbitals
                  txt='nio.txt',
                  kpts=(k, k, k),
                  xc='PBE')
e = atoms.get_potential_energy()
