#!/usr/bin/env python
"""Transport exersice

This file should do the same as pt_h2_lcao.py, but extracts the Hamiltonians
manually instead of using gpawtransport, which currently does not work
"""

from ase import Atoms
from gpaw import GPAW, Mixer, FermiDirac
from gpaw.lcao.tools import remove_pbc, get_lcao_hamiltonian, get_lead_lcao_hamiltonian
import cPickle as pickle

a = 2.41  # Pt binding lenght
b = 0.90  # H2 binding lenght
c = 1.70  # Pt-H binding lenght
L = 7.00  # width of unit cell

#####################
# Scattering region #
#####################

# Setup the Atoms for the scattering region.
atoms = Atoms('Pt5H2Pt5', pbc=(1, 0, 0), cell=[9 * a + b + 2 * c, L, L])
atoms.positions[:5, 0] = [i * a for i in range(5)]
atoms.positions[-5:, 0] = [i * a + b + 2 * c for i in range(4, 9)]
atoms.positions[5:7, 0] = [4 * a + c, 4 * a + c + b]
atoms.positions[:, 1:] = L / 2.

# Attach a GPAW calculator
calc = GPAW(h=0.3,
            xc='PBE',
            basis='szp(dzp)',
            occupations=FermiDirac(width=0.1),
            kpts=(1, 1, 1),
            mode='lcao',
            txt='pt_h2_lcao_scat.txt',
            mixer=Mixer(0.1, 5, weight=100.0),
            symmetry={'point_group': False, 'time_reversal': False})
atoms.set_calculator(calc)

atoms.get_potential_energy()  # Converge everything!
Ef = atoms.calc.get_fermi_level()

H_skMM, S_kMM = get_lcao_hamiltonian(calc)
# Only use first kpt, spin, as there are no more
H, S = H_skMM[0, 0], S_kMM[0]
H -= Ef * S
remove_pbc(atoms, H, S, 0)

# Dump the Hamiltonian and Scattering matrix to a pickle file
pickle.dump((H, S), open('scat_hs.pickle', 'wb'), 2)

########################
# Left principal layer #
########################

# Use four Pt atoms in the lead, so only take those from before
atoms = atoms[:4].copy()
atoms.set_cell([4 * a, L, L])

# Attach a GPAW calculator
calc = GPAW(h=0.3,
            xc='PBE',
            basis='szp(dzp)',
            occupations=FermiDirac(width=0.1),
            kpts=(4, 1, 1),  # More kpts needed as the x-direction is shorter
            mode='lcao',
            txt='pt_h2_lcao_llead.txt',
            mixer=Mixer(0.1, 5, weight=100.0),
            symmetry={'point_group': False, 'time_reversal': False})
atoms.set_calculator(calc)

atoms.get_potential_energy()  # Converge everything!
Ef = atoms.calc.get_fermi_level()

ibz2d_k, weight2d_k, H_skMM, S_kMM = get_lead_lcao_hamiltonian(calc)
# Only use first kpt, spin, as there are no more
H, S = H_skMM[0, 0], S_kMM[0]
H -= Ef * S

# Dump the Hamiltonian and Scattering matrix to a pickle file
pickle.dump((H, S), open('lead1_hs.pickle', 'wb'), 2)

#########################
# Right principal layer #
#########################
# This is identical to the left prinicpal layer so we don't have to do anything
# Just dump the same Hamiltonian and Scattering matrix to a pickle file
pickle.dump((H, S), open('lead2_hs.pickle', 'wb'), 2)
