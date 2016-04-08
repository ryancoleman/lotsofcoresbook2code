# This script exercsises some of the idiosyncracies
# of the descriptor class and PBLAS on a realistic
# case. See the BLACS descriptor documentation 
# in trunk/gpaw/blacs.py for some discussions of 
# these idiosyncracies.
import numpy as np

from gpaw.blacs import BlacsGrid, parallelprint
from gpaw.mpi import world, rank, size
from gpaw.utilities.scalapack import pblas_simple_gemm

gen = np.random.RandomState(42)

# simulate state-parallelization=2 and
# domain-decomposition.prod=32
B = 2
D = 32
mb = 32
grid = BlacsGrid(world, B, D)

nbands = 500
nG = 80**3

nGdesc = grid.new_descriptor(nbands, nG, nbands/B, nG/D)
nndesc = grid.new_descriptor(nbands, nbands, mb, mb)

psit_nG = gen.rand(*nGdesc.shape)
A_nn = gen.rand(*nndesc.shape)

assert nGdesc.check(psit_nG)
assert nndesc.check(A_nn)

parallelprint(world, (A_nn.shape, nndesc.shape, nndesc.lld))

pblas_simple_gemm(nGdesc, nGdesc, nndesc, psit_nG, psit_nG, A_nn,
                  transa='N', transb='T')
