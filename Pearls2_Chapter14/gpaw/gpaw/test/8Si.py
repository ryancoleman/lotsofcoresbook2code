from __future__ import print_function
from ase import Atoms
from ase.calculators.test import numeric_force
from gpaw import GPAW, FermiDirac, PoissonSolver
from gpaw.test import equal
from gpaw.xc.tools import vxc
a = 5.404
bulk = Atoms(symbols='Si8',
             positions=[(0, 0, 0.1 / a),
                        (0, 0.5, 0.5),
                        (0.5, 0, 0.5),
                        (0.5, 0.5, 0),
                        (0.25, 0.25, 0.25),
                        (0.25, 0.75, 0.75),
                        (0.75, 0.25, 0.75),
                        (0.75, 0.75, 0.25)],
             pbc=True)
bulk.set_cell((a, a, a), scale_atoms=True)
n = 20
calc = GPAW(gpts=(n, n, n),
            nbands='150%',
            occupations=FermiDirac(width=0.01),
            poissonsolver=PoissonSolver(nn='M', relax='J'),
            kpts=(2, 2, 2),
            convergence={'energy': 1e-7}
            )
bulk.set_calculator(calc)
f1 = bulk.get_forces()[0, 2]
e1 = bulk.get_potential_energy()
v_xc = vxc(calc)
print(v_xc)
niter1 = calc.get_number_of_iterations()

f2 = numeric_force(bulk, 0, 2)
print((f1, f2, f1 - f2))
equal(f1, f2, 0.005)

# Volume per atom:
vol = a**3 / 8
de = calc.get_electrostatic_corrections() / vol
print(de)
assert abs(de[0] - -2.190) < 0.001

print((e1, f1, niter1))
energy_tolerance = 0.00025
force_tolerance = 0.0001
niter_tolerance = 0
equal(e1, -46.6596470348, energy_tolerance)  # svnversion 5252
equal(f1, -1.38242356123, force_tolerance)  # svnversion 5252
