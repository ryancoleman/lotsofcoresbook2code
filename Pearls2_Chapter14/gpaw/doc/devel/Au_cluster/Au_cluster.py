#!/usr/bin/env python

from ase import Atoms
from ase.io.xyz import read_xyz

from gpaw import GPAW
from gpaw.mixer import Mixer
from gpaw import ConvergenceError
from gpaw.mpi import rank

from gpaw.eigensolvers.rmm_diis_old import RMM_DIIS

from gpaw import setup_paths

# Use setups from the $PWD and $PWD/.. first
setup_paths.insert(0, '.')
setup_paths.insert(0, '../')

atoms = read_xyz('../Au102_revised.xyz')

prefix = 'Au_cluster'
L = 32.0
atoms.set_cell((L,L,L),scale_atoms=False)
atoms.center()
atoms.set_pbc(1)
r = [1, 1, 1]
atoms = atoms.repeat(r)
n = [240 * ri for ri in r]
# nbands (>=1683) is the number of bands per cluster
nbands = 3*6*6*16 # 1728
for ri in r: nbands = nbands*ri
mixer = Mixer(beta=0.1, nmaxold=5, weight=100.0)
# the next three lines decrease memory usage
es = RMM_DIIS(keep_htpsit=False)
from gpaw.hs_operators import MatrixOperator
MatrixOperator.nblocks = 16
calc = GPAW(nbands=nbands,
            # uncomment next two lines to use lcao/sz
            #mode='lcao',
            #basis='sz',
            gpts=tuple(n),
            maxiter=5,
            width = 0.1,
            xc='LDA',
            mixer=mixer,
            eigensolver = es,
            txt=prefix + '.txt',
            )
atoms.set_calculator(calc)
from gpaw.mpi import rank
try:
    pot = atoms.get_potential_energy()
except ConvergenceError:
    pass
