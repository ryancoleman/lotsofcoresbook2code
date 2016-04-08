from ase import Atoms
from ase.optimize import QuasiNewton
from gpaw import GPAW
from gpaw.poisson import PoissonSolver

a = 6
b = a / 2

mol = Atoms('H2O',
             positions=[(b, 0.7633 + b, -0.4876 + b),
                        (b, -0.7633 + b, -0.4876 + b),
                        (b, b, 0.1219 + b)],
            cell=[a, a, a])

calc = GPAW(nbands=4,
            mode='lcao',
            basis='dzp',
            gpts=(32, 32, 32),
            poissonsolver=PoissonSolver(relax='GS', eps=1e-7))

mol.set_calculator(calc)
dyn = QuasiNewton(mol, trajectory='lcao2_h2o.traj')
dyn.run(fmax=0.05)
