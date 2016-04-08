from __future__ import print_function
import numpy as np
import ase.units as units
from gpaw import restart, GPAW
from gpaw.poisson import PoissonSolver
from gpaw.dipole_correction import DipoleCorrection

energies = []
for name in ['zero', 'periodic', 'corrected']:
    if name == 'corrected':
        calc = GPAW(name, txt=None,
                    poissonsolver=DipoleCorrection(PoissonSolver(), 2))
    else:
        calc = GPAW(name, txt=None)

    energies.append(calc.get_potential_energy())

print(energies)
assert abs(energies[1] - energies[0]) < 0.003
assert abs(energies[2] - energies[0] - 0.0409) < 0.003

efermi = calc.get_fermi_level()
calc.restore_state()
v = (calc.hamiltonian.vHt_g * units.Hartree).mean(0).mean(0)
w1 = v[0] - efermi
w2 = v[-1] - efermi
print(w1, w2)
assert abs(w1 - 4.359) < 0.01
assert abs(w2 - 2.556) < 0.01
