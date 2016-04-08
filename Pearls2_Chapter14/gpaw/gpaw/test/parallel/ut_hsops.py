from __future__ import print_function
import gc
import sys
import time
import numpy as np

try:
    # Matplotlib is not a dependency
    import matplotlib as mpl
    mpl.use('Agg')  # force the antigrain backend
except (ImportError, RuntimeError):
    mpl = None

from ase.units import Bohr
from ase.utils import devnull

from gpaw.mpi import world, distribute_cpus
from gpaw.utilities.tools import tri2full, md5_array, gram_schmidt
from gpaw.band_descriptor import BandDescriptor
from gpaw.grid_descriptor import GridDescriptor
from gpaw.kohnsham_layouts import BandLayouts
from gpaw.hs_operators import MatrixOperator
from gpaw.parameters import InputParameters
from gpaw.xc import XC
from gpaw.setup import SetupData, Setups
from gpaw.lfc import LFC

# -------------------------------------------------------------------

from gpaw.test.ut_common import ase_svnversion, shapeopt, TestCase, \
    TextTestRunner, CustomTextTestRunner, defaultTestLoader, \
    initialTestLoader, create_random_atoms, create_parsize_maxbands

memstats = False
if memstats:
    # Developer use of this feature requires ASE 3.1.0 svn.rev. 905 or later.
    assert ase_svnversion >= 905 # wasn't bug-free untill 973!
    from ase.utils.memory import MemorySingleton, MemoryStatistics

# -------------------------------------------------------------------

p = InputParameters(spinpol=False)
xc = XC(p.xc)
p.setups = dict([(symbol, SetupData(symbol, xc.name)) for symbol in 'HN'])

class UTBandParallelSetup(TestCase):
    """
    Setup a simple band parallel calculation."""

    # Number of bands
    nbands = 36

    # Spin-paired, single kpoint
    nspins = 1
    nibzkpts = 1

    # Strided or blocked groups
    parstride_bands = None

    # Mean spacing and number of grid points per axis (G x G x G)
    h = 1.0 / Bohr
    G = 20

    # Wavefunction data type
    dtype = None

    def setUp(self):
        for virtvar in ['dtype','parstride_bands']:
            assert getattr(self,virtvar) is not None, 'Virtual "%s"!' % virtvar

        parsize_domain, parsize_bands = create_parsize_maxbands(self.nbands, world.size)
        assert self.nbands % parsize_bands == 0
        comms = distribute_cpus(parsize_domain,
                                parsize_bands, self.nspins, self.nibzkpts)
        domain_comm, kpt_comm, band_comm, block_comm = \
            [comms[name] for name in 'dkbK']
        self.block_comm = block_comm

        # Set up band descriptor:
        self.bd = BandDescriptor(self.nbands, band_comm, self.parstride_bands)

        # Set up grid descriptor:
        res, ngpts = shapeopt(100, self.G**3, 3, 0.2)
        cell_c = self.h * np.array(ngpts)
        pbc_c = (True, False, True)
        self.gd = GridDescriptor(ngpts, cell_c, pbc_c, domain_comm, parsize_domain)

        # Create Kohn-Sham layouts for these band and grid descriptors:
        self.ksl = self.create_kohn_sham_layouts()

        # What to do about kpoints?
        self.kpt_comm = kpt_comm

    def tearDown(self):
        del self.bd, self.gd, self.ksl, self.kpt_comm, self.block_comm

    def create_kohn_sham_layouts(self):
        return BandLayouts(self.gd, self.bd, self.block_comm, self.dtype)

    # =================================

    def verify_comm_sizes(self):
        if world.size == 1:
            return
        comm_sizes = tuple([comm.size for comm in [world, self.bd.comm, \
                                                   self.gd.comm, self.kpt_comm]])
        self._parinfo =  '%d world, %d band, %d domain, %d kpt' % comm_sizes
        self.assertEqual((self.nspins*self.nibzkpts) % self.kpt_comm.size, 0)

    def verify_kohn_sham_layouts(self):
        # TODO do more here :)
        self.assertFalse(self.ksl.using_blacs)
        self.assertTrue(self.ksl.bd is self.bd)
        self.assertTrue(self.ksl.gd is self.gd)

    def verify_band_stride_related(self):
        # Verify that (q1+q2)%B-q1 falls in ]-B;Q[ where Q=B//2+1
        for B in range(1,256):
            Q = B//2+1
            #dqs = []
            #for q1 in range(B):
            #    for q2 in range(Q):
            #        dq = (q1+q2)%B-q1
            #        dqs.append(dq)
            #dqs = np.array(dqs)
            q1s = np.arange(B)[:,np.newaxis]
            q2s = np.arange(Q)[np.newaxis,:]
            dqs = (q1s+q2s)%B-q1s
            self.assertEqual(dqs.min(), -B+1)
            self.assertEqual(dqs.max(), Q-1)

    def verify_band_indexing_consistency(self):
        for n in range(self.bd.nbands):
            band_rank, myn = self.bd.who_has(n)
            self.assertEqual(self.bd.global_index(myn, band_rank), n)

        for band_rank in range(self.bd.comm.size):
            for myn in range(self.bd.mynbands):
                n = self.bd.global_index(myn, band_rank)
                self.assertTrue(self.bd.who_has(n) == (band_rank, myn))

    def verify_band_ranking_consistency(self):
        rank_n = self.bd.get_band_ranks()

        for band_rank in range(self.bd.comm.size):
            my_band_indices = self.bd.get_band_indices(band_rank)
            matches = np.argwhere(rank_n == band_rank).ravel()
            self.assertTrue((matches == my_band_indices).all())

            for myn in range(self.bd.mynbands):
                n = self.bd.global_index(myn, band_rank)
                self.assertEqual(my_band_indices[myn], n)


class UTBandParallelSetup_Blocked(UTBandParallelSetup):
    __doc__ = UTBandParallelSetup.__doc__
    parstride_bands = False
    dtype = float

class UTBandParallelSetup_Strided(UTBandParallelSetup):
    __doc__ = UTBandParallelSetup.__doc__
    parstride_bands = True
    dtype = float
# -------------------------------------------------------------------

def record_memory(wait=0.1):
    assert gc.collect() == 0, 'Uncollected garbage!'
    world.barrier()
    time.sleep(wait)
    mem = MemoryStatistics()
    time.sleep(wait)
    world.barrier()
    return mem

def create_memory_info(mem1, mem2, vmkey='VmData'):
    dm = np.array([(mem2-mem1)[vmkey]], dtype=float)
    dm_r = np.empty(world.size, dtype=float)
    world.all_gather(dm, dm_r)
    return dm_r, ','.join(map(lambda v: '%8.4f MB' % v, dm_r/1024**2.))

# -------------------------------------------------------------------

class UTConstantWavefunctionSetup(UTBandParallelSetup):
    __doc__ = UTBandParallelSetup.__doc__ + """
    The pseudo wavefunctions are constants normalized to their band index."""

    allocated = False
    blocking = None
    async = None
    
    def setUp(self):
        UTBandParallelSetup.setUp(self)
        for virtvar in ['dtype','blocking','async']:
            assert getattr(self,virtvar) is not None, 'Virtual "%s"!' % virtvar

        # Create randomized atoms
        self.atoms = create_random_atoms(self.gd)

        # Do we agree on the atomic positions?
        pos_ac = self.atoms.get_positions()
        pos_rac = np.empty((world.size,)+pos_ac.shape, pos_ac.dtype)
        world.all_gather(pos_ac, pos_rac)
        if (pos_rac-pos_rac[world.rank,...][np.newaxis,...]).any():
            raise RuntimeError('Discrepancy in atomic positions detected.')

        # Create setups for atoms
        self.Z_a = self.atoms.get_atomic_numbers()
        self.setups = Setups(self.Z_a, p.setups, p.basis,
                             p.lmax, xc)

        # Create atomic projector overlaps
        spos_ac = self.atoms.get_scaled_positions() % 1.0
        self.rank_a = self.gd.get_ranks_from_positions(spos_ac)
        self.pt = LFC(self.gd, [setup.pt_j for setup in self.setups],
                      dtype=self.dtype)
        self.pt.set_positions(spos_ac)

        if memstats:
            # Hack to scramble heap usage into steady-state level
            HEAPSIZE = 25 * 1024**2
            for i in range(100):
                data = np.empty(np.random.uniform(0, HEAPSIZE // 8), float)
                del data
            self.mem_pre = record_memory()
            self.mem_alloc = None
            self.mem_test = None

        # Stuff for pseudo wave functions and projections
        if self.dtype == complex:
            self.gamma = 1j**(1.0/self.nbands)
        else:
            self.gamma = 1.0

        self.psit_nG = None
        self.P_ani = None
        self.Qeff_a = None
        self.Qtotal = None

        self.allocate()

    def tearDown(self):
        UTBandParallelSetup.tearDown(self)
        del self.P_ani, self.psit_nG
        del self.pt, self.setups, self.atoms
        if memstats:
            self.print_memory_summary()
            del self.mem_pre, self.mem_alloc, self.mem_test
        self.allocated = False

    def print_memory_summary(self):
        if not memstats:
            raise RuntimeError('No memory statistics were recorded!')

        if world.rank == 0:
            sys.stdout.write('\n')
            sys.stdout.flush()

        dm_r, dminfo = create_memory_info(MemorySingleton(), self.mem_pre)
        if world.rank == 0:
            print('overhead: %s -> %8.4f MB' % (dminfo, dm_r.sum()/1024**2.))

        dm_r, dminfo = create_memory_info(self.mem_pre, self.mem_alloc)
        if world.rank == 0:
            print('allocate: %s -> %8.4f MB' % (dminfo, dm_r.sum()/1024**2.))

        dm_r, dminfo = create_memory_info(self.mem_alloc, self.mem_test)
        if world.rank == 0:
            print('test-use: %s -> %8.4f MB' % (dminfo, dm_r.sum()/1024**2.))

    def allocate(self):
        """
        Allocate constant wavefunctions and their projections according to::

                          /          \     _____    i*phase*(n-m)
           <psi |psi > = ( 1 + Q      ) * V m*n  * e
               n    n'    \     total/

        """
        if self.allocated:
            raise RuntimeError('Already allocated!')

        self.allocate_wavefunctions()
        self.allocate_projections()

        # XXX DEBUG disables projection contributions
        if False:
            for a,P_ni in self.P_ani.items():
                P_ni[:] = 0.
                self.Qeff_a[a] = 0.
            self.Qtotal = 0.
            self.Z_a[:] = 2**63-1
        # XXX DEBUG

        self.Qtotal = np.empty(1, dtype=float)
        self.Qtotal[:] = np.sum([Qeff for Qeff in self.Qeff_a.values()])
        self.gd.comm.sum(self.Qtotal)

        band_indices = np.arange(self.nbands).astype(self.dtype)
        z = self.gamma**band_indices * band_indices**0.5
        self.S0_nn = (1. + self.Qtotal) * np.outer(z.conj(), z)

        self.allocated = True

        if memstats:
            self.mem_alloc = record_memory()

    def allocate_wavefunctions(self):
        """
        Allocate constant pseudo wavefunctions according to::

             ~    ~        _____    i*phase*(n-m)
           <psi |psi > =  V m*n  * e
               n    n'

        """
        if self.allocated:
            raise RuntimeError('Already allocated!')

        # Fill in wave functions
        gpts_c = self.gd.get_size_of_global_array()
        self.psit_nG = self.gd.empty(self.bd.mynbands, self.dtype)
        for myn, psit_G in enumerate(self.psit_nG):
            n = self.bd.global_index(myn)
            # Fill psit_nG: | psit_n > = exp(i*phase*n) * sqrt(n) / sqrt(V)
            psit_G[:] = self.gamma**n * n**0.5 / (self.gd.dv * gpts_c.prod())**0.5

    def allocate_projections(self):
        """
        Construct dummy projection of pseudo wavefunction according to::

           ___
           \     ~   ~a    a   ~a  ~             1     _____    i*phase*(n-m)
            )  <psi |p > dO   <p |psi > =  +/-  --- * V m*n  * e
           /___    n  i    ii'  i'   n'          Z
            ii'                                   a

        """
        if self.allocated:
            raise RuntimeError('Already allocated!')

        # Fill in projector overlaps
        my_band_indices = self.bd.get_band_indices()
        my_atom_indices = np.argwhere(self.gd.comm.rank == self.rank_a).ravel()

        # Holm-Nielsen check:
        natoms = len(self.atoms)
        assert (self.gd.comm.sum(float(sum(my_atom_indices))) ==
                natoms * (natoms - 1) // 2)

        # Check that LFC agrees with us:
        self.assertEqual(len(my_atom_indices), len(self.pt.my_atom_indices))
        for a1, a2 in zip(my_atom_indices, self.pt.my_atom_indices):
            self.assertEqual(a1, a2)

        self.Qeff_a = {}
        self.P_ani = self.pt.dict(self.bd.mynbands)
        for a in my_atom_indices:
            ni = self.setups[a].ni
            # Fill P_ni: <p_i | psit_n > = beta_i * exp(i*phase*n) * sqrt(n)
            #
            #  |  ____                   |
            #  |  \        *    a        |      1
            #  |   )   beta   dO   beta  |  =  ----
            #  |  /___     i    ij     j |      Z
            #  |    ij                   |       a
            #
            # Substitution by linear transformation: beta_i ->  dO_ij alpha_j,
            # where we start out with some initial non-constant vector:
            alpha_i = np.exp(-np.arange(ni).astype(self.dtype)/ni)
            try:
                # Try Cholesky decomposition dO_ii = L_ii * L_ii^dag
                L_ii = np.linalg.cholesky(self.setups[a].dO_ii)
                alpha_i /= np.vdot(alpha_i, alpha_i)**0.5
                beta_i = np.linalg.solve(L_ii.T.conj(), alpha_i)
            except np.linalg.LinAlgError:
                # Eigenvector decomposition dO_ii = V_ii * W_ii * V_ii^dag
                W_i, V_ii = np.linalg.eigh(self.setups[a].dO_ii)
                alpha_i /= np.abs(np.vdot(alpha_i, 
                                          np.dot(np.diag(W_i), alpha_i)))**0.5
                beta_i = np.linalg.solve(V_ii.T.conj(), alpha_i)

            # Normalize according to plus/minus charge
            beta_i /= self.Z_a[a]**0.5
            self.Qeff_a[a] = np.vdot(beta_i, np.dot(self.setups[a].dO_ii, \
                                                    beta_i)).real
            self.P_ani[a][:] = np.outer(self.gamma**my_band_indices \
                                        * my_band_indices**0.5, beta_i)

    def check_and_plot(self, A_nn, A0_nn, digits, keywords=''):
        # Construct fingerprint of input matrices for comparison
        fingerprint = np.array([md5_array(A_nn, numeric=True),
                                md5_array(A0_nn, numeric=True)])

        # Compare fingerprints across all processors
        fingerprints = np.empty((world.size, 2), np.int64)
        world.all_gather(fingerprint, fingerprints)
        if fingerprints.ptp(0).any():
            raise RuntimeError('Distributed matrices are not identical!')

        # If assertion fails, catch temporarily while plotting, then re-raise
        try:
            self.assertAlmostEqual(np.abs(A_nn-A0_nn).max(), 0, digits)
        except AssertionError:
            if world.rank == 0 and mpl is not None:
                from matplotlib.figure import Figure
                fig = Figure()
                ax = fig.add_axes([0.0, 0.1, 1.0, 0.83])
                ax.set_title(self.__class__.__name__)
                im = ax.imshow(np.abs(A_nn-A0_nn), interpolation='nearest')
                fig.colorbar(im)
                fig.text(0.5, 0.05, 'Keywords: ' + keywords, \
                    horizontalalignment='center', verticalalignment='top')

                from matplotlib.backends.backend_agg import FigureCanvasAgg
                img = 'ut_hsops_%s_%s.png' % (self.__class__.__name__, \
                    '_'.join(keywords.split(',')))
                FigureCanvasAgg(fig).print_figure(img.lower(), dpi=90)
            raise

    def get_optimal_number_of_blocks(self, blocking='fast'):
        """Estimate the optimal number of blocks for band parallelization.

        The number of blocks determines how many parallel send/receive 
        operations are performed, as well as the added memory footprint 
        of the required send/receive buffers.

        ``blocking``  ``nblocks``      Description
        ============  =============    ========================================
        'fast'        ``1``            Heavy on memory, more accurate and fast.
        'light'       ``mynbands``     Light on memory, less accurate and slow.
        'intdiv'      ``...``          First integer divisible value 
        'nonintdiv'   ``...``          Some non-integer divisible cases
        """

        #if self.bd.comm.size == 1:
        #    return 1

        if blocking == 'fast':
            return 1
        elif blocking == 'light':
            return self.bd.mynbands
        elif blocking == 'intdiv':
            # Find first value of nblocks that leads to integer
            # divisible mybands / nblock. This is very like to be 
            # 2 but coded here for the general case
            nblocks = 2
            while self.bd.mynbands % nblocks != 0:
                nblocks +=1
            return nblocks
        elif blocking == 'nonintdiv1':
            # Find first value of nblocks that leads to non-integer
            # divisible mynbands / nblock that is less than M
            nblocks = 2
            M = self.bd.mynbands // nblocks
            while self.bd.mynbands % nblocks < M:
                nblocks += 1
                M = self.bd.mynbands // nblocks
            return nblocks
        elif blocking == 'nonintdiv2':
            # Find first value of nblocks that leads to non-integer
            # divisible mynbands / nblock that is less than M
            nblocks = 2
            M = self.bd.mynbands // nblocks
            while self.bd.mynbands % nblocks > M:
                nblocks += 1
                M = self.mynbands // nblocks
            return nblocks
        else:
            nblocks = blocking
            assert self.bd.mynbands // nblocks > 0
            return nblocks

    # =================================

    def test_contents_wavefunction(self):
        # Integrate diagonal brakets of pseudo wavefunctions
        gpts_c = self.gd.get_size_of_global_array()

        intpsit_myn = self.bd.empty(dtype=self.dtype)
        for myn, psit_G in enumerate(self.psit_nG):
            n = self.bd.global_index(myn)
            intpsit_myn[myn] = np.vdot(psit_G, psit_G) * self.gd.dv
        self.gd.comm.sum(intpsit_myn)

        if memstats:
            self.mem_test = record_memory()

        my_band_indices = self.bd.get_band_indices()
        self.assertAlmostEqual(np.abs(intpsit_myn-my_band_indices).max(), 0, 9)

        intpsit_n = self.bd.collect(intpsit_myn, broadcast=True)
        self.assertAlmostEqual(np.abs(intpsit_n-np.arange(self.nbands)).max(), 0, 9)

    def test_contents_projection(self):
        # Distribute inverse effective charges to everybody in domain
        all_Qeff_a = np.empty(len(self.atoms), dtype=float)
        for a,rank in enumerate(self.rank_a):
            if rank == self.gd.comm.rank:
                Qeff = np.array([self.Qeff_a[a]])
            else:
                Qeff = np.empty(1, dtype=float)
            self.gd.comm.broadcast(Qeff, rank)
            all_Qeff_a[a] = Qeff

        # Check absolute values consistency of inverse effective charges
        self.assertAlmostEqual(np.abs(1./self.Z_a-np.abs(all_Qeff_a)).max(), 0, 9)

        # Check sum of inverse effective charges against total
        self.assertAlmostEqual(all_Qeff_a.sum(), self.Qtotal, 9)

        # Make sure that we all agree on inverse effective charges
        fingerprint = np.array([md5_array(all_Qeff_a, numeric=True)])
        all_fingerprints = np.empty(world.size, fingerprint.dtype)
        world.all_gather(fingerprint, all_fingerprints)
        if all_fingerprints.ptp(0).any():
            raise RuntimeError('Distributed eff. charges are not identical!')

    def test_overlaps_hermitian(self):
        # Set up Hermitian overlap operator:
        S = lambda x: x
        dS = lambda a, P_ni: np.dot(P_ni, self.setups[a].dO_ii)
        nblocks = self.get_optimal_number_of_blocks(self.blocking)
        overlap = MatrixOperator(self.ksl, nblocks, self.async, True)
        S_nn = overlap.calculate_matrix_elements(self.psit_nG, \
            self.P_ani, S, dS).T.copy() # transpose to get <psit_m|A|psit_n>
        tri2full(S_nn, 'U') # upper to lower...

        if self.bd.comm.rank == 0:
            self.gd.comm.broadcast(S_nn, 0)
        self.bd.comm.broadcast(S_nn, 0)

        if memstats:
            self.mem_test = record_memory()

        self.check_and_plot(S_nn, self.S0_nn, 9, 'overlaps,hermitian')

    def test_overlaps_nonhermitian(self):
        alpha = np.random.normal(size=1).astype(self.dtype)
        if self.dtype == complex:
            alpha += 1j*np.random.normal(size=1)
        world.broadcast(alpha, 0)

        # Set up non-Hermitian overlap operator:
        S = lambda x: alpha*x
        dS = lambda a, P_ni: np.dot(alpha*P_ni, self.setups[a].dO_ii)
        nblocks = self.get_optimal_number_of_blocks(self.blocking)
        overlap = MatrixOperator(self.ksl, nblocks, self.async, False)
        S_nn = overlap.calculate_matrix_elements(self.psit_nG, \
            self.P_ani, S, dS).T.copy() # transpose to get <psit_m|A|psit_n>

        if self.bd.comm.rank == 0:
            self.gd.comm.broadcast(S_nn, 0)
        self.bd.comm.broadcast(S_nn, 0)

        if memstats:
            self.mem_test = record_memory()

        self.check_and_plot(S_nn, alpha*self.S0_nn, 9, 'overlaps,nonhermitian')

    def test_trivial_cholesky(self):
        # Known starting point of SI_nn = <psit_m|S+alpha*I|psit_n>
        I_nn = np.eye(*self.S0_nn.shape)
        alpha = 1e-3 # shift eigenvalues away from zero
        SI_nn = self.S0_nn + alpha * I_nn

        # Try Cholesky decomposition SI_nn = L_nn * L_nn^dag
        L_nn = np.linalg.cholesky(SI_nn)
        # |psit_n> -> C_nn |psit_n> , C_nn^(-1) = L_nn^dag
        # <psit_m|SI|psit_n> -> <psit_m|C_nn^dag SI C_nn|psit_n> = diag(W_n)
        C_nn = np.linalg.inv(L_nn.T.conj())

        # Set up Hermitian overlap operator:
        S = lambda x: x
        dS = lambda a, P_ni: np.dot(P_ni, self.setups[a].dO_ii)
        nblocks = self.get_optimal_number_of_blocks(self.blocking)
        overlap = MatrixOperator(self.ksl, nblocks, self.async, True)
        self.psit_nG = overlap.matrix_multiply(C_nn.T.copy(), self.psit_nG, self.P_ani)
        D_nn = overlap.calculate_matrix_elements(self.psit_nG, \
            self.P_ani, S, dS).T.copy() # transpose to get <psit_m|A|psit_n>
        tri2full(D_nn, 'U') # upper to lower...

        if self.bd.comm.rank == 0:
            self.gd.comm.broadcast(D_nn, 0)
        self.bd.comm.broadcast(D_nn, 0)

        if memstats:
            self.mem_test = record_memory()

        # D_nn = C_nn^dag * S_nn * C_nn = I_nn - alpha * C_nn^dag * C_nn
        D0_nn = I_nn - alpha * np.dot(C_nn.T.conj(), C_nn)
        self.check_and_plot(D_nn, D0_nn, 6, 'trivial,cholesky') #XXX precision

    def test_trivial_diagonalize(self):
        # Known starting point of S_nn = <psit_m|S|psit_n>
        S_nn = self.S0_nn

        # Eigenvector decomposition S_nn = V_nn * W_nn * V_nn^dag
        # Utilize the fact that they are analytically known (cf. Maple)
        band_indices = np.arange(self.nbands)
        V_nn = np.eye(self.nbands).astype(self.dtype)
        if self.dtype == complex:
            V_nn[1:,1] = np.conj(self.gamma)**band_indices[1:] * band_indices[1:]**0.5
            V_nn[1,2:] = -self.gamma**band_indices[1:-1] * band_indices[2:]**0.5
        else:
            V_nn[2:,1] = band_indices[2:]**0.5
            V_nn[1,2:] = -band_indices[2:]**0.5

        W_n = np.zeros(self.nbands).astype(self.dtype)
        W_n[1] = (1. + self.Qtotal) * self.nbands * (self.nbands - 1) / 2.

        # Find the inverse basis
        Vinv_nn = np.linalg.inv(V_nn)

        # Test analytical eigenvectors for consistency against analytical S_nn
        D_nn = np.dot(Vinv_nn, np.dot(S_nn, V_nn))
        self.assertAlmostEqual(np.abs(D_nn.diagonal()-W_n).max(), 0, 8)
        self.assertAlmostEqual(np.abs(np.tril(D_nn, -1)).max(), 0, 4)
        self.assertAlmostEqual(np.abs(np.triu(D_nn, 1)).max(), 0, 4)
        del Vinv_nn, D_nn

        # Perform Gram Schmidt orthonormalization for diagonalization
        # |psit_n> -> C_nn |psit_n>, using orthonormalized basis Q_nn
        # <psit_m|S|psit_n> -> <psit_m|C_nn^dag S C_nn|psit_n> = diag(W_n)
        # using S_nn = V_nn * W_nn * V_nn^(-1) = Q_nn * W_nn * Q_nn^dag
        C_nn = V_nn.copy()
        gram_schmidt(C_nn)
        self.assertAlmostEqual(np.abs(np.dot(C_nn.T.conj(), C_nn) \
                                      - np.eye(self.nbands)).max(), 0, 6)

        # Set up Hermitian overlap operator:
        S = lambda x: x
        dS = lambda a, P_ni: np.dot(P_ni, self.setups[a].dO_ii)
        nblocks = self.get_optimal_number_of_blocks(self.blocking)
        overlap = MatrixOperator(self.ksl, nblocks, self.async, True)
        self.psit_nG = overlap.matrix_multiply(C_nn.T.copy(), self.psit_nG, self.P_ani)
        D_nn = overlap.calculate_matrix_elements(self.psit_nG, \
            self.P_ani, S, dS).T.copy() # transpose to get <psit_m|A|psit_n>
        tri2full(D_nn, 'U') # upper to lower...

        if self.bd.comm.rank == 0:
            self.gd.comm.broadcast(D_nn, 0)
        self.bd.comm.broadcast(D_nn, 0)

        if memstats:
            self.mem_test = record_memory()

        # D_nn = C_nn^dag * S_nn * C_nn = W_n since Q_nn^dag = Q_nn^(-1)
        D0_nn = np.dot(C_nn.T.conj(), np.dot(S_nn, C_nn))
        self.assertAlmostEqual(np.abs(D0_nn-np.diag(W_n)).max(), 0, 9)
        self.check_and_plot(D_nn, D0_nn, 9, 'trivial,diagonalize')

    def test_multiply_randomized(self):
        # Known starting point of S_nn = <psit_m|S|psit_n>
        S_nn = self.S0_nn

        if self.dtype == complex:
            C_nn = np.random.uniform(size=self.nbands**2) * \
                np.exp(1j*np.random.uniform(0,2*np.pi,size=self.nbands**2))
        else:
            C_nn = np.random.normal(size=self.nbands**2)
        C_nn = C_nn.reshape((self.nbands,self.nbands)) / np.linalg.norm(C_nn,2)
        world.broadcast(C_nn, 0)

        # Set up Hermitian overlap operator:
        S = lambda x: x
        dS = lambda a, P_ni: np.dot(P_ni, self.setups[a].dO_ii)
        nblocks = self.get_optimal_number_of_blocks(self.blocking)
        overlap = MatrixOperator(self.ksl, nblocks, self.async, True)
        self.psit_nG = overlap.matrix_multiply(C_nn.T.copy(), self.psit_nG, self.P_ani)
        D_nn = overlap.calculate_matrix_elements(self.psit_nG, \
            self.P_ani, S, dS).T.copy() # transpose to get <psit_m|A|psit_n>
        tri2full(D_nn, 'U') # upper to lower...

        if self.bd.comm.rank == 0:
            self.gd.comm.broadcast(D_nn, 0)
        self.bd.comm.broadcast(D_nn, 0)

        if memstats:
            self.mem_test = record_memory()

        # D_nn = C_nn^dag * S_nn * C_nn
        D0_nn = np.dot(C_nn.T.conj(), np.dot(S_nn, C_nn))
        self.check_and_plot(D_nn, D0_nn, 9, 'multiply,randomized')

    def test_multiply_nonhermitian(self):
        alpha = np.random.normal(size=1).astype(self.dtype)
        if self.dtype == complex:
            alpha += 1j*np.random.normal(size=1)
        world.broadcast(alpha, 0)

        # Known starting point of S_nn = <psit_m|S|psit_n>
        S_nn = alpha*self.S0_nn

        if self.dtype == complex:
            C_nn = np.random.uniform(size=self.nbands**2) * \
                np.exp(1j*np.random.uniform(0,2*np.pi,size=self.nbands**2))
        else:
            C_nn = np.random.normal(size=self.nbands**2)
        C_nn = C_nn.reshape((self.nbands,self.nbands)) / np.linalg.norm(C_nn,2)
        world.broadcast(C_nn, 0)

        # Set up non-Hermitian overlap operator:
        S = lambda x: alpha*x
        dS = lambda a, P_ni: np.dot(alpha*P_ni, self.setups[a].dO_ii)
        nblocks = self.get_optimal_number_of_blocks(self.blocking)
        overlap = MatrixOperator(self.ksl, nblocks, self.async, False)
        self.psit_nG = overlap.matrix_multiply(C_nn.T.copy(), self.psit_nG, self.P_ani)
        D_nn = overlap.calculate_matrix_elements(self.psit_nG, \
            self.P_ani, S, dS).T.copy() # transpose to get <psit_m|A|psit_n>

        if self.bd.comm.rank == 0:
            self.gd.comm.broadcast(D_nn, 0)
        self.bd.comm.broadcast(D_nn, 0)

        if memstats:
            self.mem_test = record_memory()

        # D_nn = C_nn^dag * S_nn * C_nn
        D0_nn = np.dot(C_nn.T.conj(), np.dot(S_nn, C_nn))
        self.check_and_plot(D_nn, D0_nn, 9, 'multiply,nonhermitian')


# -------------------------------------------------------------------

def UTConstantWavefunctionFactory(dtype, parstride_bands, blocking, async):
    sep = '_'
    classname = 'UTConstantWavefunctionSetup' \
    + sep + {float:'Float', complex:'Complex'}[dtype] \
    + sep + {False:'Blocked', True:'Strided'}[parstride_bands] \
    + sep + {'fast':'Fast', 'light':'Light', 
             'intdiv':'Intdiv', 'nonintdiv1':'Nonintdiv1',
             'nonintdiv2':'Nonintdiv2'}[blocking] \
    + sep + {False:'Synchronous', True:'Asynchronous'}[async]
    class MetaPrototype(UTConstantWavefunctionSetup, object):
        __doc__ = UTConstantWavefunctionSetup.__doc__
        dtype = dtype
        parstride_bands = parstride_bands
        blocking = blocking
        async = async
    MetaPrototype.__name__ = classname
    return MetaPrototype

# -------------------------------------------------------------------

if __name__ in ['__main__', '__builtin__']:
    # We may have been imported by test.py, if so we should redirect to logfile
    if __name__ == '__builtin__':
        testrunner = CustomTextTestRunner('ut_hsops.log', verbosity=2)
    else:
        stream = (world.rank == 0) and sys.stdout or devnull
        testrunner = TextTestRunner(stream=stream, verbosity=2)

    parinfo = []
    # Initial Verification only tests case : dtype = float
    for test in [UTBandParallelSetup_Blocked, UTBandParallelSetup_Strided]:
        info = ['', test.__name__, test.__doc__.strip('\n'), '']
        testsuite = initialTestLoader.loadTestsFromTestCase(test)
        map(testrunner.stream.writeln, info)
        testresult = testrunner.run(testsuite)
        assert testresult.wasSuccessful(), 'Initial verification failed!'
        parinfo.extend(['    Parallelization options: %s' % tci._parinfo for \
                        tci in testsuite._tests if hasattr(tci, '_parinfo')])
    parinfo = np.unique(np.sort(parinfo)).tolist()

    testcases = []
    for dtype in [float, complex]:
        for parstride_bands in [False, True]:
            for blocking in ['fast', 'light', 'intdiv',   
                             'nonintdiv1', 'nonintdiv2']: 
                for async in [False, True]:
                    testcases.append(UTConstantWavefunctionFactory(dtype, \
                        parstride_bands, blocking, async))

    for test in testcases:
        info = ['', test.__name__, test.__doc__.strip('\n')] + parinfo + ['']
        testsuite = defaultTestLoader.loadTestsFromTestCase(test)
        map(testrunner.stream.writeln, info)
        testresult = testrunner.run(testsuite)
        # Provide feedback on failed tests if imported by test.py
        if __name__ == '__builtin__' and not testresult.wasSuccessful():
            raise SystemExit('Test failed. Check ut_hsops.log for details.')

