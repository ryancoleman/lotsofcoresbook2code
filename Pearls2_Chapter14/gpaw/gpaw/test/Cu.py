from __future__ import print_function
import os
from ase import Atoms
from ase.units import Hartree
from gpaw import GPAW
from gpaw.test import equal, gen
import gpaw.mpi as mpi

# Generate non-scalar-relativistic setup for Cu:
gen('Cu', scalarrel=False)

a = 8.0
c = a / 2
Cu = Atoms('Cu', [(c, c, c)], magmoms=[1],
           cell=(a, a, a), pbc=0)

calc = GPAW(h=0.2, lmax=0)# basis='sz')
Cu.set_calculator(calc)
e = Cu.get_potential_energy()
niter = calc.get_number_of_iterations()

e_4s_major = calc.get_eigenvalues(spin=0)[5] / Hartree
e_3d_minor = calc.get_eigenvalues(spin=1)[4] / Hartree
print(mpi.rank, e_4s_major, e_3d_minor)

#
# The reference values are from:
#
#   http://physics.nist.gov/PhysRefData/DFTdata/Tables/29Cu.html
#
if mpi.rank == 0:
    print(e_4s_major - e_3d_minor, -0.184013 - -0.197109)
    assert abs(e_4s_major - e_3d_minor - (-0.184013 - -0.197109)) < 0.001

    print(e, niter)
    energy_tolerance = 0.0005
    niter_tolerance = 0
    equal(e, -0.271504, energy_tolerance)
