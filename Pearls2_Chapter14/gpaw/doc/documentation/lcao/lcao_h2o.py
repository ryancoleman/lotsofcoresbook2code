from ase import Atoms
from ase.optimize import QuasiNewton
from gpaw import GPAW

a = 6
b = a / 2

mol = Atoms('H2O',
            [(b, 0.7633 + b, -0.4876 + b),
             (b, -0.7633 + b, -0.4876 + b),
             (b, b, 0.1219 + b)],
            cell=[a, a, a])

calc = GPAW(nbands=4,
            h=0.2,
            mode='lcao',
            basis='dzp')

mol.set_calculator(calc)
dyn = QuasiNewton(mol, trajectory='lcao_h2o.traj')
dyn.run(fmax=0.05)
