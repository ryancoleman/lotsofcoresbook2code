from __future__ import print_function
from gpaw import GPAW, FermiDirac
from gpaw.mpi import world, size, rank
from gpaw.lrtddft2 import LrTDDFTindexed, lr_communicators
from gpaw.test import equal
from ase.atoms import Atoms
from glob import glob
import os
import numpy as np

debug = False
use_hdf5 = True

try:
    import _gpaw_hdf5
    restart_file = 'Al2_gs.hdf5'
except ImportError:
    restart_file = 'Al2_gs.gpw'

d = 2.563
atoms = Atoms('Al2', positions=((0, 0, 0),
                                (0, 0, d))
             )
atoms.center(4.0)
calc = GPAW(h=0.24, eigensolver='cg', basis='dzp',
            occupations=FermiDirac(width=0.01),
            convergence={'eigenstates': 4.0e-5, 'density' : 1.0e-2, 
                         'bands' : 'all'},
            nbands=20)
atoms.set_calculator(calc)
atoms.get_potential_energy()
calc.write(restart_file, mode='all')

# Try to run parallel over eh-pairs
if size % 2 == 0:
    eh_size = 2
    domain_size = size // eh_size
else:
    eh_size = 1
    domain_size = size

dd_comm, eh_comm = lr_communicators(world, domain_size, eh_size)

calc = GPAW(restart_file,
            communicator=dd_comm)
de = 3.0
lr = LrTDDFTindexed('Al2_lri',
                    calc=calc,
                    xc = 'PBE',
                    max_energy_diff=de,
                    eh_communicator=eh_comm
                   )
lr.calculate_excitations()
e0_1 = np.sqrt(lr.evalues[0])
e1_1 = np.sqrt(lr.evalues[-1])
# Continue with larger de
de = 4.5
lr = LrTDDFTindexed('Al2_lri',
                    calc=calc,
                    xc = 'PBE',
                    max_energy_diff=de,
                    eh_communicator=eh_comm
                   )
lr.calculate_excitations()
e0_2 = np.sqrt(lr.evalues[0])
e1_2 = np.sqrt(lr.evalues[-1])
# Continue with smaller de
de = 2.5
lr = LrTDDFTindexed('Al2_lri',
                    calc=calc,
                    xc = 'PBE',
                    max_energy_diff=de,
                    eh_communicator=eh_comm
                   )
lr.calculate_excitations()
e0_3 = np.sqrt(lr.evalues[0])
e1_3 = np.sqrt(lr.evalues[-1])

if debug and rank == 0:
    print(e0_1, e1_1)
    print(e0_2, e1_2)
    print(e0_3, e1_3)

tol = 1.0e-8
equal(e0_1, 0.00105074187176, tol)
equal(e0_1, e0_2, tol)
equal(e0_1, e0_3, tol)
equal(e1_1, 0.183188157301, tol)
equal(e1_2, 0.194973135812, tol)
equal(e1_3, 0.120681529342, tol)

# Remove the unused output files
if rank == 0:
    files = glob('*.ready_rows.*')
    files += glob('*.LR_info')
    files += glob('*.log.*')
    files += glob('*.K_matrix.*')
    files += glob('*.KS_singles')
    for file in files:
        os.remove(file)
