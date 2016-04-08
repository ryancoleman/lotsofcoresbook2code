""" test run for DFTB+ calculator 
    the tolerance is extremely loose, beause different sk files 
    give different results

"""

from ase.test.dftb import installed
assert installed()

from ase import Atoms
from ase.visualize import view
from ase.calculators.dftb import Dftb
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.lattice.surface import *
from ase.lattice.surface import add_adsorbate
import os, sys

h = 1.85
d = 1.10

slab = fcc111('Ni', size=(2,2,3), vacuum=10.0)
calc1=\
Dftb(label='slab',
kpts=[2,2,1],
Hamiltonian_MaxAngularMomentum_ = '',
Hamiltonian_MaxAngularMomentum_Ni = '"d"',
Hamiltonian_SCC='YES'
)
slab.set_calculator(calc1)
dyn = QuasiNewton(slab, trajectory='slab.traj')
dyn.run(fmax=0.05)
e_slab = slab.get_potential_energy()
os.system('rm dftb_in.hsd')
molecule = Atoms('2N', positions=[(0., 0., 0.), (0., 0., d)])
calc2=\
Dftb(label='n2',
Hamiltonian_MaxAngularMomentum_ = '',
Hamiltonian_MaxAngularMomentum_N = '"p"',
Hamiltonian_SCC='YES'
)
molecule.set_calculator(calc2)
dyn = QuasiNewton(molecule, trajectory='n2.traj')
dyn.run(fmax=0.05)
e_N2 = molecule.get_potential_energy()

slab2 = slab
add_adsorbate(slab2, molecule, h, 'ontop')
constraint = FixAtoms(mask=[a.symbol != 'N' for a in slab2])
slab2.set_constraint(constraint)
calc3=\
Dftb(label='slab2',
kpts=[2,2,1],
Hamiltonian_MaxAngularMomentum_ = '',
Hamiltonian_MaxAngularMomentum_N = '"p"',
Hamiltonian_MaxAngularMomentum_Ni = '"d"',
Hamiltonian_SCC='YES'
)
slab2.set_calculator(calc3)
dyn = QuasiNewton(slab2, trajectory='slab2.traj')
dyn.run(fmax=0.05)

adsorption_energy = e_slab + e_N2 - slab2.get_potential_energy()
assert abs(adsorption_energy + 0.4227415) < 0.3
files = ['band.out', 'detailed.out', 'dftb_in.hsd', 'dftb_pin.hsd', \
    'geo_end.gen', 'geo_end.xyz', 'results.tag', \
    'n2.traj', 'n2.traj.bak', 'n2.out', \
    'slab.traj', 'slab.traj.bak', 'slab.out', \
    'slab2.traj', 'slab2.traj.bak', 'slab2.out' \
]

for file in files:
    try:
        os.remove(file)
    except OSError:
        pass

