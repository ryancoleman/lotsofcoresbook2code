import numpy as np
from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from gpaw.test import equal

a = 5.475
calc = GPAW(h=0.24,
            kpts=(4, 4, 4),
            occupations=FermiDirac(width=0.0),
            nbands=5)
atoms = bulk('Si', 'diamond', a=a)
atoms.set_calculator(calc)
E = atoms.get_potential_energy()
equal(E, -11.8092858, 0.0002)
niter = calc.get_number_of_iterations()

equal(atoms.calc.get_fermi_level(), 5.17751284, 0.005)
homo, lumo = calc.get_homo_lumo()
equal(lumo - homo, 1.11445025, 0.002)

calc.write('si_primitive.gpw', 'all')
calc = GPAW('si_primitive.gpw',
            parallel={'domain': 1, 'band': 1},
            idiotproof=False,
            txt=None)
from gpaw.xc.hybridk import HybridXC
pbe0 = HybridXC('PBE0', alpha=5.0)
calc.get_xc_difference(pbe0)
