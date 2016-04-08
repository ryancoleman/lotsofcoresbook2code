from __future__ import print_function
from time import time
import sys
import numpy as np
from gpaw import parsize_domain, parsize_bands
from gpaw.band_descriptor import BandDescriptor
from gpaw.grid_descriptor import GridDescriptor
from gpaw.kohnsham_layouts import BandLayouts
from gpaw.mpi import world, distribute_cpus
from gpaw.utilities import gcd
from gpaw.utilities.lapack import inverse_cholesky
from gpaw.hs_operators import MatrixOperator

G = 120  # number of grid points (G x G x G)
N = 2000  # number of bands
repeats = 20

try:
    N = int(sys.argv[1])
    K = int(sys.argv[2])
except (IndexError, ValueError):
    N = 6
    K = 3
    repeats = 3

# B: number of band groups
# D: number of domains
if parsize_bands is None:
    if parsize_domain is None:
        B = gcd(N, world.size)
        D = world.size // B
    else:
        B = world.size // np.prod(parsize_domain)
        D = parsize_domain
else:
    B = parsize_bands
    D = world.size // B

M = N // B     # number of bands per group
assert M * B == N, 'M=%d, B=%d, N=%d' % (M,B,N)

h = 0.2        # grid spacing
a = h * G      # side length of box
assert np.prod(D) * B == world.size, 'D=%s, B=%d, W=%d' % (D,B,world.size)

# Set up communicators:
comms = distribute_cpus(parsize_domain=D,
                        parsize_bands=B,
                        nspins=1, nibzkpts=1)
domain_comm, kpt_comm, band_comm, block_comm = \
    [comms[name] for name in ['d', 'k', 'b', 'K']]


assert kpt_comm.size == 1
if world.rank == 0:
    print('MPI: %d domains, %d band groups' % (domain_comm.size, band_comm.size))

# Set up band and grid descriptors:
bd = BandDescriptor(N, band_comm, False)
gd = GridDescriptor((G, G, G), (a, a, a), True, domain_comm, parsize=D)
ksl = BandLayouts(gd, bd, block_comm, float)

# Random wave functions:
psit_mG = gd.empty(M)
for m in range(M):
    np.random.seed(world.rank * M + m)
    psit_mG[m] = np.random.uniform(-0.5, 0.5, tuple(gd.n_c))
if world.rank == 0:
    print('Size of wave function array:', psit_mG.shape)
P_ani = {0: psit_mG[:, :2, 0, 0].copy(),
         1: psit_mG[:, -1, -1, -3:].copy()}
if 0:
    X = M // K
    assert K * X == M
    if G**3 // D // K * K < G**3 // D:
        X += 1
    print(X)
    work1_xG = gd.empty(X)
    work2_xG = gd.empty(X)

def run(psit_mG):
    overlap = MatrixOperator(ksl, K)
    if 0:
        overlap.work1_xG = work1_xG
        overlap.work2_xG = work2_xG
    #S_nn = np.empty((N, N))
    def S(x):
        return x
    dS_aii = {0: np.ones((2, 2)) * 0.123, 1: np.ones((3, 3)) * 0.321}
    def dS(a, P_ni):
        return np.dot(P_ni, dS_aii[a])
    S_nn = overlap.calculate_matrix_elements(psit_mG, P_ani, S, dS)

    t1 = time()
    if world.rank == 0:
        print(S_nn.round(5))
        inverse_cholesky(S_nn)
    C_nn = S_nn
    t2 = time()

    if world.rank == 0:
        print('Cholesky Time %f' % (t2-t1))
        
    # Distribute matrix:
    world.broadcast(C_nn, 0)

    psit_mG = overlap.matrix_multiply(C_nn, psit_mG, P_ani)

    if world.rank == 0:
        print('Made it past matrix multiply')

    # Check:
    S_nn = overlap.calculate_matrix_elements(psit_mG, P_ani, S, dS)

    assert not(P_ani[0] - psit_mG[:, :2, 0, 0]).round(10).any()
    assert not(P_ani[1] - psit_mG[:, -1, -1, -3:]).round(10).any()

    if world.rank == 0:
        for n in range(N):
            assert abs(S_nn[n, n] - 1.0) < 1e-10
            assert not S_nn[n + 1:, n].round(10).any()

    return psit_mG

ta = time()

# Do twenty iterations
for x in range(repeats):
    psit_mG = run(psit_mG)

tb = time()

if world.rank == 0:
    print('Total Time %f' % (tb -ta))
    
