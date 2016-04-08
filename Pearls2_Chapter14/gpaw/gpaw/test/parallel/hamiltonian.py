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
from gpaw.fd_operators import Laplace

G = 120  # number of grid points (G x G x G)
N = 2000  # number of bands
repeats = 1

try:
    N = int(sys.argv[1])
    J = int(sys.argv[2])
except (IndexError, ValueError):
    N = 6
    J = 3
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

kin = Laplace(gd, -0.5, 2).apply
vt_G = gd.empty()
vt_G.fill(0.567)

def run(psit_mG):
    overlap = MatrixOperator(ksl, J)
    def H(psit_xG):
        Htpsit_xG = np.empty_like(psit_xG)
        kin(psit_xG, Htpsit_xG)
        for psit_G, y_G in zip(psit_xG, Htpsit_xG):
            y_G += vt_G * psit_G
        return Htpsit_xG
    dH_aii = {0: np.ones((2, 2)) * 0.123, 1: np.ones((3, 3)) * 0.321}
    def dH(a, P_ni):
        return np.dot(P_ni, dH_aii[a])
    H_nn = overlap.calculate_matrix_elements(psit_mG, P_ani, H, dH)

    t1 = time()
    if world.rank == 0:
        eps_n, H_nn = np.linalg.eigh(H_nn)
        H_nn = np.ascontiguousarray(H_nn.T)
    t2 = time()

    if world.rank == 0:
        print('Diagonalization Time %f' % (t2-t1))
        print(eps_n)

    # Distribute matrix:
    world.broadcast(H_nn, 0)

    psit_mG = overlap.matrix_multiply(H_nn, psit_mG, P_ani)

    if world.rank == 0:
        print('Made it past matrix multiply')

    # Check:
    assert not(P_ani[0] - psit_mG[:, :2, 0, 0]).round(10).any()
    assert not(P_ani[1] - psit_mG[:, -1, -1, -3:]).round(10).any()

    H_nn = overlap.calculate_matrix_elements(psit_mG, P_ani, H, dH)
    if world.rank == 0:
        for n in range(N):
            assert abs(H_nn[n, n] - eps_n[n]) < 1.5e-8
            assert not H_nn[n + 1:, n].round(8).any()

    return psit_mG

ta = time()

for x in range(repeats):
    psit_mG = run(psit_mG)

tb = time()

if world.rank == 0:
    print('Total Time %f' % (tb -ta))
