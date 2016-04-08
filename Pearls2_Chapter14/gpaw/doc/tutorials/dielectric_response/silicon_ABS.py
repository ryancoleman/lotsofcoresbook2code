# Refer to G. Kresse, Phys. Rev. B 73, 045112 (2006)
# for comparison of macroscopic and microscopic dielectric constant
# and absorption peaks.
from __future__ import print_function

from ase.lattice import bulk
from ase.parallel import paropen
from gpaw import GPAW, FermiDirac
from gpaw.response.df import DielectricFunction

# Ground state calculation
a = 5.431
atoms = bulk('Si', 'diamond', a=a)

calc = GPAW(mode='pw',
            kpts={'density': 5.0, 'gamma': True},
            xc='LDA',
            occupations=FermiDirac(0.001))  # Use small FD smearing

atoms.set_calculator(calc)
atoms.get_potential_energy()  # Get ground state density

# Restart Calculation with fixed density and dense kpoint sampling
calc.set(kpts={'density': 15.0, 'gamma': False},  # Dense kpoint sampling
         fixdensity=True)
atoms.get_potential_energy()
calc.diagonalize_full_hamiltonian(nbands=70)  # Diagonalize Hamiltonian
calc.write('si_large.gpw', 'all')  # Write wavefunctions

# Getting absorption spectrum
df = DielectricFunction(calc='si_large.gpw',
                        eta=0.05,
                        domega0=0.02,
                        ecut=150)
df.get_dielectric_function(filename='si_abs.csv')

# Getting macroscopic constant
df = DielectricFunction(calc='si_large.gpw',
                        frequencies=[0.0],
                        hilbert=False,
                        eta=0.0001,
                        ecut=150,
                        )

epsNLF, epsLF = df.get_macroscopic_dielectric_constant()

# Make table
epsrefNLF = 14.08  # From [1] in top
epsrefLF = 12.66  # From [1] in top

f = paropen('mac_eps.csv', 'w')
print(' , Without LFE, With LFE', file=f)
print('%s, %.6f, %.6f' % ('GPAW-linear response', epsNLF, epsLF), file=f)
print('%s, %.6f, %.6f' % ('[1]', epsrefNLF, epsrefLF), file=f)
print('%s, %.6f, %.6f' % ('Exp.', 11.90, 11.90), file=f)
f.close()
