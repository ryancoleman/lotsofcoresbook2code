
import sys
import numpy as np

from ase.utils import devnull

from gpaw import debug
from gpaw.mpi import world
from gpaw.utilities.tools import tri2full
from gpaw.hs_operators import MatrixOperator
from gpaw.utilities import compiled_with_sl
from gpaw.utilities.scalapack import scalapack_set
from gpaw.blacs import Redistributor
from gpaw.kohnsham_layouts import BlacsBandLayouts 
if 0:  # causes numpy doctests failures - exact formatting is expected!
    np.set_printoptions(linewidth=168) #XXX large xterm width

# -------------------------------------------------------------------

from gpaw.test.ut_common import ase_svnversion, TextTestRunner, \
    CustomTextTestRunner, defaultTestLoader, initialTestLoader

memstats = False
if memstats:
    # Developer use of this feature requires ASE 3.1.0 svn.rev. 905 or later.
    assert ase_svnversion >= 905 # wasn't bug-free untill 973!
    from ase.utils.memory import MemorySingleton, MemoryStatistics

from gpaw.test.parallel.ut_hsops import UTBandParallelSetup, \
                                        UTConstantWavefunctionSetup

# -------------------------------------------------------------------

class UTBandParallelBlacsSetup(UTBandParallelSetup):
    """
    Setup a simple band parallel calculation using BLACS."""

    def create_kohn_sham_layouts(self):
        # Find BLACS parameters for Kohn-Sham layouts
        cpus = self.bd.comm.size * self.gd.comm.size
        self.mcpus = int(cpus**0.5)
        self.ncpus = cpus//self.mcpus
        return BlacsBandLayouts(self.gd, self.bd, self.block_comm, self.dtype, self.mcpus, self.ncpus, 6)

    # =================================

    def verify_comm_sizes(self):
        if world.size == 1:
            return
        comm_sizes = tuple([comm.size for comm in [world, self.bd.comm, \
                                                   self.gd.comm, self.kpt_comm]])
        comm_sizes += (self.mcpus, self.ncpus)
        self._parinfo =  '%d world, %d band, %d domain, %d kpt, %d x %d BLACS' % comm_sizes
        self.assertEqual((self.nspins*self.nibzkpts) % self.kpt_comm.size, 0)

    def verify_kohn_sham_layouts(self):
        # TODO do more here :)
        self.assertTrue(self.ksl.using_blacs)
        self.assertTrue(self.ksl.bd is self.bd)
        self.assertTrue(self.ksl.gd is self.gd)


class UTBandParallelBlacsSetup_Blocked(UTBandParallelBlacsSetup):
    __doc__ = UTBandParallelBlacsSetup.__doc__
    parstride_bands = False
    dtype = float

class UTBandParallelBlacsSetup_Strided(UTBandParallelSetup):
    __doc__ = UTBandParallelBlacsSetup.__doc__
    parstride_bands = True
    dtype = float

# -------------------------------------------------------------------

class UTConstantWavefunctionBlacsSetup(UTConstantWavefunctionSetup,
                                       UTBandParallelBlacsSetup):
    __doc__ = UTBandParallelBlacsSetup.__doc__ + """
    The pseudo wavefunctions are constants normalized to their band index."""

    # =================================

    def test_overlaps_hermitian(self):
        # Set up Hermitian overlap operator:
        S = lambda x: x
        dS = lambda a, P_ni: np.dot(P_ni, self.setups[a].dO_ii)
        nblocks = self.get_optimal_number_of_blocks(self.blocking)
        overlap = MatrixOperator(self.ksl, nblocks, self.async, True)
        S_nn = overlap.calculate_matrix_elements(self.psit_nG, self.P_ani, S, dS)

        if memstats:
            self.mem_test = record_memory()

        S_NN = self.ksl.nndescriptor.collect_on_master(S_nn)
        if self.bd.comm.rank == 0 and self.gd.comm.rank == 0:
            assert S_NN.shape == (self.bd.nbands,) * 2
            S_NN = S_NN.T.copy() # Fortran -> C indexing
            tri2full(S_NN, 'U') # upper to lower...
        else:
            assert S_NN.nbytes == 0
            S_NN = np.empty((self.bd.nbands,) * 2, dtype=S_NN.dtype)

        if self.bd.comm.rank == 0:
            self.gd.comm.broadcast(S_NN, 0)
        self.bd.comm.broadcast(S_NN, 0)

        self.check_and_plot(S_NN, self.S0_nn, 9, 'overlaps,hermitian')

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
        if 0: #XXX non-hermitian case so Nn2nn not just uplo='L' but rather 'G'
            blockcomm = self.ksl.nndescriptor.blacsgrid.comm
            self.ksl.Nn2nn = Redistributor(blockcomm, self.ksl.Nndescriptor,
                                           self.ksl.nndescriptor)
        S_nn = overlap.calculate_matrix_elements(self.psit_nG, self.P_ani, S, dS)

        if memstats:
            self.mem_test = record_memory()

        S_NN = self.ksl.nndescriptor.collect_on_master(S_nn)
        if self.bd.comm.rank == 0 and self.gd.comm.rank == 0:
            assert S_NN.shape == (self.bd.nbands,) * 2
            S_NN = S_NN.T.copy() # Fortran -> C indexing
        else:
            assert S_NN.nbytes == 0
            S_NN = np.empty((self.bd.nbands,) * 2, dtype=S_NN.dtype)

        if self.bd.comm.rank == 0:
            self.gd.comm.broadcast(S_NN, 0)
        self.bd.comm.broadcast(S_NN, 0)

        self.check_and_plot(S_NN, alpha*self.S0_nn, 9, 'overlaps,nonhermitian')

    def test_trivial_cholesky(self):
        # Set up Hermitian overlap operator:
        S = lambda x: x
        dS = lambda a, P_ni: np.dot(P_ni, self.setups[a].dO_ii)
        nblocks = self.get_optimal_number_of_blocks(self.blocking)
        overlap = MatrixOperator(self.ksl, nblocks, self.async, True)
        S_nn = overlap.calculate_matrix_elements(self.psit_nG, self.P_ani, S, dS)

        # Known starting point of SI_nn = <psit_m|S+alpha*I|psit_n>
        I_nn = self.ksl.nndescriptor.empty(dtype=S_nn.dtype)
        scalapack_set(self.ksl.nndescriptor, I_nn, 0.0, 1.0, 'L')
        alpha = 1e-3 # shift eigenvalues away from zero

        C_nn = S_nn + alpha * I_nn
        self.ksl.nndescriptor.inverse_cholesky(C_nn, 'L')
        self.psit_nG = overlap.matrix_multiply(C_nn, self.psit_nG, self.P_ani)
        D_nn = overlap.calculate_matrix_elements(self.psit_nG, self.P_ani, S, dS)

        D_NN = self.ksl.nndescriptor.collect_on_master(D_nn)
        if self.bd.comm.rank == 0 and self.gd.comm.rank == 0:
            assert D_NN.shape == (self.bd.nbands,) * 2
            D_NN = D_NN.T.copy() # Fortran -> C indexing
            tri2full(D_NN, 'U') # upper to lower..
        else:
            assert D_NN.nbytes == 0
            D_NN = np.empty((self.bd.nbands,) * 2, dtype=D_NN.dtype)

        if self.bd.comm.rank == 0:
            self.gd.comm.broadcast(D_NN, 0)
        self.bd.comm.broadcast(D_NN, 0)

        # D_NN = C_NN^dag * S_NN * C_NN = I_NN - alpha * C_NN^dag * C_NN
        I_NN = np.eye(self.bd.nbands)
        C0_NN = np.linalg.inv(np.linalg.cholesky(self.S0_nn + alpha*I_NN).T.conj())
        D0_NN = I_NN - alpha * np.dot(C0_NN.T.conj(), C0_NN)
        self.check_and_plot(D_NN, D0_NN, 6, 'trivial,cholesky') #XXX precision

    def test_trivial_diagonalize(self): #XXX XXX XXX
        # Known starting point of S_nn = <psit_m|S|psit_n>
        S_nn = self.S0_nn

        # Eigenvector decomposition S_nn = V_nn * W_nn * V_nn^dag
        # Utilize the fact that they are analytically known (cf. Maple)
        W_n = np.zeros(self.nbands).astype(self.dtype)
        W_n[1] = (1. + self.Qtotal) * self.nbands * (self.nbands - 1) / 2.

        # Set up Hermitian overlap operator:
        S = lambda x: x
        dS = lambda a, P_ni: np.dot(P_ni, self.setups[a].dO_ii)
        nblocks = self.get_optimal_number_of_blocks(self.blocking)
        overlap = MatrixOperator(self.ksl, nblocks, self.async, True)
        S_nn = overlap.calculate_matrix_elements(self.psit_nG, self.P_ani, S, dS)

        eps_N = self.bd.empty(global_array=True) # XXX dtype?
        C_nn = self.ksl.nndescriptor.empty(dtype=S_nn.dtype)
        self.ksl.nndescriptor.diagonalize_dc(S_nn, C_nn, eps_N, 'L')
        self.assertAlmostEqual(np.abs(np.sort(eps_N)-np.sort(W_n)).max(), 0, 9)

        #eps_n = self.bd.empty()
        #self.bd.distribute(eps_N, eps_n) # XXX only blocked groups, right?

        # Rotate wavefunctions to diagonalize the overlap
        self.psit_nG = overlap.matrix_multiply(C_nn, self.psit_nG, self.P_ani)

        # Recaulculate the overlap matrix, which should now be diagonal
        D_nn = overlap.calculate_matrix_elements(self.psit_nG, self.P_ani, S, dS)

        D_NN = self.ksl.nndescriptor.collect_on_master(D_nn)
        if self.bd.comm.rank == 0 and self.gd.comm.rank == 0:
            assert D_NN.shape == (self.bd.nbands,) * 2
            D_NN = D_NN.T.copy() # Fortran -> C indexing
            tri2full(D_NN, 'U') # upper to lower...
        else:
            assert D_NN.nbytes == 0
            D_NN = np.empty((self.bd.nbands,) * 2, dtype=D_NN.dtype)

        if self.bd.comm.rank == 0:
            self.gd.comm.broadcast(D_NN, 0)
        self.bd.comm.broadcast(D_NN, 0)

        D0_NN = np.diag(eps_N)
        self.check_and_plot(D_NN, D0_NN, 9, 'trivial,diagonalize')

    def test_multiply_randomized(self):
        # Known starting point of S_nn = <psit_m|S|psit_n>
        S_NN = self.S0_nn

        if self.dtype == complex:
            C_NN = np.random.uniform(size=self.nbands**2) * \
                np.exp(1j*np.random.uniform(0,2*np.pi,size=self.nbands**2))
        else:
            C_NN = np.random.normal(size=self.nbands**2)
        C_NN = C_NN.reshape((self.nbands,self.nbands)) / np.linalg.norm(C_NN,2)
        world.broadcast(C_NN, 0)

        # Set up Hermitian overlap operator:
        S = lambda x: x
        dS = lambda a, P_ni: np.dot(P_ni, self.setups[a].dO_ii)
        nblocks = self.get_optimal_number_of_blocks(self.blocking)
        overlap = MatrixOperator(self.ksl, nblocks, self.async, True)

        if self.bd.comm.rank == 0 and self.gd.comm.rank == 0:
            assert C_NN.shape == (self.bd.nbands,) * 2
            tmp_NN = C_NN.T.copy() # C -> Fortran indexing
        else:
            tmp_NN = self.ksl.nndescriptor.as_serial().empty(dtype=C_NN.dtype)
        C_nn = self.ksl.nndescriptor.distribute_from_master(tmp_NN)

        self.psit_nG = overlap.matrix_multiply(C_nn, self.psit_nG, self.P_ani)
        D_nn = overlap.calculate_matrix_elements(self.psit_nG, self.P_ani, S, dS)

        if memstats:
            self.mem_test = record_memory()

        D_NN = self.ksl.nndescriptor.collect_on_master(D_nn)
        if self.bd.comm.rank == 0 and self.gd.comm.rank == 0:
            assert D_NN.shape == (self.bd.nbands,) * 2
            D_NN = D_NN.T.copy() # Fortran -> C indexing
            tri2full(D_NN, 'U') # upper to lower...
        else:
            assert D_NN.nbytes == 0
            D_NN = np.empty((self.bd.nbands,) * 2, dtype=D_NN.dtype)

        if self.bd.comm.rank == 0:
            self.gd.comm.broadcast(D_NN, 0)
        self.bd.comm.broadcast(D_NN, 0)

        # D_nn = C_nn^dag * S_nn * C_nn
        D0_NN = np.dot(C_NN.T.conj(), np.dot(S_NN, C_NN))
        self.check_and_plot(D_NN, D0_NN, 9, 'multiply,randomized')

    def test_multiply_nonhermitian(self):
        alpha = np.random.normal(size=1).astype(self.dtype)
        if self.dtype == complex:
            alpha += 1j*np.random.normal(size=1)
        world.broadcast(alpha, 0)

        # Known starting point of S_nn = <psit_m|S|psit_n>
        S_NN = alpha*self.S0_nn

        if self.dtype == complex:
            C_NN = np.random.uniform(size=self.nbands**2) * \
                np.exp(1j*np.random.uniform(0,2*np.pi,size=self.nbands**2))
        else:
            C_NN = np.random.normal(size=self.nbands**2)
        C_NN = C_NN.reshape((self.nbands,self.nbands)) / np.linalg.norm(C_NN,2)
        world.broadcast(C_NN, 0)

        # Set up Hermitian overlap operator:
        S = lambda x: alpha*x
        dS = lambda a, P_ni: np.dot(alpha*P_ni, self.setups[a].dO_ii)
        nblocks = self.get_optimal_number_of_blocks(self.blocking)
        overlap = MatrixOperator(self.ksl, nblocks, self.async, False)

        if 0: #XXX non-hermitian case so Nn2nn not just uplo='L' but rather 'G'
            blockcomm = self.ksl.nndescriptor.blacsgrid.comm
            self.ksl.Nn2nn = Redistributor(blockcomm, self.ksl.Nndescriptor,
                                           self.ksl.nndescriptor)

        if self.bd.comm.rank == 0 and self.gd.comm.rank == 0:
            assert C_NN.shape == (self.bd.nbands,) * 2
            tmp_NN = C_NN.T.copy() # C -> Fortran indexing
        else:
            tmp_NN = self.ksl.nndescriptor.as_serial().empty(dtype=C_NN.dtype)
        C_nn = self.ksl.nndescriptor.distribute_from_master(tmp_NN)

        self.psit_nG = overlap.matrix_multiply(C_nn, self.psit_nG, self.P_ani)
        D_nn = overlap.calculate_matrix_elements(self.psit_nG, self.P_ani, S, dS)

        if memstats:
            self.mem_test = record_memory()

        D_NN = self.ksl.nndescriptor.collect_on_master(D_nn)
        if self.bd.comm.rank == 0 and self.gd.comm.rank == 0:
            assert D_NN.shape == (self.bd.nbands,) * 2
            D_NN = D_NN.T.copy() # Fortran -> C indexing
        else:
            assert D_NN.nbytes == 0
            D_NN = np.empty((self.bd.nbands,) * 2, dtype=D_NN.dtype)

        if self.bd.comm.rank == 0:
            self.gd.comm.broadcast(D_NN, 0)
        self.bd.comm.broadcast(D_NN, 0)

        # D_nn = C_nn^dag * S_nn * C_nn
        D0_NN = np.dot(C_NN.T.conj(), np.dot(S_NN, C_NN))
        self.check_and_plot(D_NN, D0_NN, 9, 'multiply,nonhermitian')


# -------------------------------------------------------------------

def UTConstantWavefunctionFactory(dtype, parstride_bands, blocking, async):
    sep = '_'
    classname = 'UTConstantWavefunctionBlacsSetup' \
    + sep + {float:'Float', complex:'Complex'}[dtype] \
    + sep + {False:'Blocked', True:'Strided'}[parstride_bands] \
    + sep + {'fast':'Fast', 'light':'Light', 
             'intdiv':'Intdiv', 'nonintdiv1': 'Nonintdiv1', 
             'nonintdiv2':'Nonintdiv2'}[blocking] \
    + sep + {False:'Synchronous', True:'Asynchronous'}[async]
    class MetaPrototype(UTConstantWavefunctionBlacsSetup, object):
        __doc__ = UTConstantWavefunctionBlacsSetup.__doc__
        dtype = dtype
        parstride_bands = parstride_bands
        blocking = blocking
        async = async
    MetaPrototype.__name__ = classname
    return MetaPrototype

# -------------------------------------------------------------------

if __name__ in ['__main__', '__builtin__'] and compiled_with_sl():
    # We may have been imported by test.py, if so we should redirect to logfile
    if __name__ == '__builtin__':
        testrunner = CustomTextTestRunner('ut_hsblacs.log', verbosity=2)
    else:
        stream = (world.rank == 0) and sys.stdout or devnull
        testrunner = TextTestRunner(stream=stream, verbosity=2)

    parinfo = []
    # Initial Verification only tests case : dtype = float
    for test in [UTBandParallelBlacsSetup_Blocked]: #, UTBandParallelBlacsSetup_Strided]:
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
        for parstride_bands in [False]: #XXX [False, True]:
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
            raise SystemExit('Test failed. Check ut_hsblacs.log for details.')

