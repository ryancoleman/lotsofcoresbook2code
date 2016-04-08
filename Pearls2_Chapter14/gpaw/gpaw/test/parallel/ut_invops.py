
import sys
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
from gpaw.utilities.tools import md5_array
from gpaw.utilities.gauss import gaussian_wave
from gpaw.band_descriptor import BandDescriptor
from gpaw.grid_descriptor import GridDescriptor
from gpaw.kpt_descriptor import KPointDescriptor
from gpaw.kohnsham_layouts import BandLayouts
from gpaw.parameters import InputParameters
from gpaw.xc import XC
from gpaw.setup import SetupData, Setups
from gpaw.wavefunctions.base import WaveFunctions
from gpaw.wavefunctions.fd import FDWaveFunctions
from gpaw.fd_operators import Laplace # required but not really used
from gpaw.pair_overlap import GridPairOverlap, ProjectorPairOverlap

# -------------------------------------------------------------------

from gpaw.test.ut_common import ase_svnversion, shapeopt, TestCase, \
    TextTestRunner, CustomTextTestRunner, defaultTestLoader, \
    initialTestLoader, create_random_atoms, create_parsize_minbands

# -------------------------------------------------------------------

p = InputParameters(spinpol=False)
xc = XC(p.xc)
p.setups = dict([(symbol, SetupData(symbol, xc.name)) for symbol in 'HO'])

class UTDomainParallelSetup(TestCase):
    """
    Setup a simple domain parallel calculation."""

    # Number of bands
    nbands = 1

    # Spin-paired, single kpoint
    nspins = 1
    nibzkpts = 1

    # Mean spacing and number of grid points per axis (G x G x G)
    h = 0.25 / Bohr
    G = 48

    # Type of boundary conditions employed
    boundaries = None

    def setUp(self):
        for virtvar in ['boundaries']:
            assert getattr(self,virtvar) is not None, 'Virtual "%s"!' % virtvar

        parsize_domain, parsize_bands = create_parsize_minbands(self.nbands, world.size)
        assert self.nbands % np.prod(parsize_bands) == 0
        comms = distribute_cpus(parsize_domain,
                                parsize_bands, self.nspins, self.nibzkpts)
        domain_comm, kpt_comm, band_comm, block_comm = \
            [comms[name] for name in 'dkbK']

        self.block_comm = block_comm

        # Set up band descriptor:
        self.bd = BandDescriptor(self.nbands, band_comm)

        # Set up grid descriptor:
        res, ngpts = shapeopt(300, self.G**3, 3, 0.2)
        cell_c = self.h * np.array(ngpts)
        pbc_c = {'zero'    : False, \
                 'periodic': True, \
                 'mixed'   : (True, False, True)}[self.boundaries]
        self.gd = GridDescriptor(ngpts, cell_c, pbc_c, domain_comm, parsize_domain)

        # What to do about kpoints?
        self.kpt_comm = kpt_comm

    def tearDown(self):
        del self.bd, self.gd, self.kpt_comm, self.block_comm

    # =================================

    def verify_comm_sizes(self):
        if world.size == 1:
            return
        comm_sizes = tuple([comm.size for comm in [world, self.bd.comm, \
                                                   self.gd.comm, self.kpt_comm]])
        self._parinfo =  '%d world, %d band, %d domain, %d kpt' % comm_sizes
        self.assertEqual(self.nbands % self.bd.comm.size, 0)
        self.assertEqual((self.nspins*self.nibzkpts) % self.kpt_comm.size, 0)


class UTDomainParallelSetup_Zero(UTDomainParallelSetup):
    __doc__ = UTDomainParallelSetup.__doc__
    boundaries = 'zero'

class UTDomainParallelSetup_Periodic(UTDomainParallelSetup):
    __doc__ = UTDomainParallelSetup.__doc__
    boundaries = 'periodic'

class UTDomainParallelSetup_Mixed(UTDomainParallelSetup):
    __doc__ = UTDomainParallelSetup.__doc__
    boundaries = 'mixed'

# -------------------------------------------------------------------

# Helper functions/classes here

class FDWFS(FDWaveFunctions):
    
    def __init__(self, gd, bd, kd, setups, block_comm, dtype): # override constructor

        assert kd.comm.size == 1

        WaveFunctions.__init__(self, gd, 1, setups, bd, dtype, world,
                               kd, None)
        self.kin = Laplace(gd, -0.5, dtype=dtype)
        self.diagksl = None
        self.orthoksl = BandLayouts(gd, bd, block_comm, dtype)
        self.initksl = None
        self.overlap = None
        self.rank_a = None

    def allocate_arrays_for_projections(self, my_atom_indices): # no alloc
        pass

    def collect_projections(self, P_ani):
        if self.gd.comm.size == 1 and self.bd.comm.size == 1:
            return np.concatenate([P_ni.T for P_ni in P_ani.values()])

        assert len(self.kpt_u) == 1
        self.kpt_u[0].P_ani = P_ani
        all_P_ni = WaveFunctions.collect_projections(self, 0, 0)
        if self.world.rank == 0:
            P_In = all_P_ni.T.copy()
        else:
            nproj = sum([setup.ni for setup in self.setups])
            P_In = np.empty((nproj, self.bd.nbands), self.pt.dtype)
        self.world.broadcast(P_In, 0)
        return P_In
        
# -------------------------------------------------------------------

class UTGaussianWavefunctionSetup(UTDomainParallelSetup):
    __doc__ = UTDomainParallelSetup.__doc__ + """
    The pseudo wavefunctions are moving gaussians centered around each atom."""

    allocated = False
    dtype = None

    # Default arguments for scaled Gaussian wave
    _sigma0 = 2.0 #0.75
    _k0_c = 2*np.pi*np.array([1/5., 1/3., 0.])

    def setUp(self):
        UTDomainParallelSetup.setUp(self)

        for virtvar in ['dtype']:
            assert getattr(self,virtvar) is not None, 'Virtual "%s"!' % virtvar

        # Create randomized atoms
        self.atoms = create_random_atoms(self.gd, 5) # also tested: 10xNH3/BDA

        # XXX DEBUG START
        if False:
            from ase import view
            view(self.atoms*(1+2*self.gd.pbc_c))
        # XXX DEBUG END

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

        # K-point descriptor
        bzk_kc = np.array([[0, 0, 0]], dtype=float)
        self.kd = KPointDescriptor(bzk_kc, 1)
        self.kd.set_symmetry(self.atoms, self.setups)
        self.kd.set_communicator(self.kpt_comm)
        
        # Create gamma-point dummy wavefunctions
        self.wfs = FDWFS(self.gd, self.bd, self.kd, self.setups,
                         self.block_comm, self.dtype)
        
        spos_ac = self.atoms.get_scaled_positions() % 1.0
        self.wfs.set_positions(spos_ac)
        self.pt = self.wfs.pt # XXX shortcut

        ## Also create pseudo partial waveves
        #from gpaw.lfc import LFC
        #self.phit = LFC(self.gd, [setup.phit_j for setup in self.setups], \
        #                self.kpt_comm, dtype=self.dtype)
        #self.phit.set_positions(spos_ac)

        self.r_cG = None
        self.buf_G = None
        self.psit_nG = None

        self.allocate()

    def tearDown(self):
        UTDomainParallelSetup.tearDown(self)
        del self.r_cG, self.buf_G, self.psit_nG
        del self.pt, self.setups, self.atoms
        self.allocated = False

    def allocate(self):
        self.r_cG = self.gd.get_grid_point_coordinates()
        cell_cv = self.atoms.get_cell() / Bohr
        assert np.abs(cell_cv-self.gd.cell_cv).max() < 1e-9
        center_c = 0.5*cell_cv.diagonal()

        self.buf_G = self.gd.empty(dtype=self.dtype)
        self.psit_nG = self.gd.empty(self.bd.mynbands, dtype=self.dtype)

        for myn,psit_G in enumerate(self.psit_nG):
            n = self.bd.global_index(myn)
            psit_G[:] = self.get_scaled_gaussian_wave(center_c, scale=10+2j*n)

            k_c = 2*np.pi*np.array([1/2., -1/7., 0.])
            for pos_c in self.atoms.get_positions() / Bohr:
                sigma = self._sigma0/(1+np.sum(pos_c**2))**0.5
                psit_G += self.get_scaled_gaussian_wave(pos_c, sigma, k_c, n+5j)

        self.allocated = True

    def get_scaled_gaussian_wave(self, pos_c, sigma=None, k_c=None, scale=None):
        if sigma is None:
            sigma = self._sigma0

        if k_c is None:
            k_c = self._k0_c

        if scale is None:
            A = None
        else:
            # 4*pi*int(exp(-r^2/(2*w^2))^2*r^2, r=0...infinity)= w^3*pi^(3/2)
            # = scale/A^2 -> A = scale*(sqrt(Pi)*w)^(-3/2) hence int -> scale^2
            A = scale/(sigma*(np.pi)**0.5)**1.5

        return gaussian_wave(self.r_cG, pos_c, sigma, k_c, A, self.dtype, self.buf_G)

    def check_and_plot(self, P_ani, P0_ani, digits, keywords=''):
        # Collapse into viewable matrices
        P_In = self.wfs.collect_projections(P_ani)
        P0_In = self.wfs.collect_projections(P0_ani)

        # Construct fingerprint of input matrices for comparison
        fingerprint = np.array([md5_array(P_In, numeric=True),
                                md5_array(P0_In, numeric=True)])

        # Compare fingerprints across all processors
        fingerprints = np.empty((world.size, 2), np.int64)
        world.all_gather(fingerprint, fingerprints)
        if fingerprints.ptp(0).any():
            raise RuntimeError('Distributed matrices are not identical!')

        # If assertion fails, catch temporarily while plotting, then re-raise
        try:
            self.assertAlmostEqual(np.abs(P_In-P0_In).max(), 0, digits)
        except AssertionError:
            if world.rank == 0 and mpl is not None:
                from matplotlib.figure import Figure
                fig = Figure()
                ax = fig.add_axes([0.0, 0.1, 1.0, 0.83])
                ax.set_title(self.__class__.__name__)
                im = ax.imshow(np.abs(P_In-P0_In), interpolation='nearest')
                fig.colorbar(im)
                fig.text(0.5, 0.05, 'Keywords: ' + keywords, \
                    horizontalalignment='center', verticalalignment='top')

                from matplotlib.backends.backend_agg import FigureCanvasAgg
                img = 'ut_invops_%s_%s.png' % (self.__class__.__name__, \
                    '_'.join(keywords.split(',')))
                FigureCanvasAgg(fig).print_figure(img.lower(), dpi=90)
            raise

    # =================================

    def test_projection_linearity(self):
        kpt = self.wfs.kpt_u[0]
        Q_ani = self.pt.dict(self.bd.mynbands)
        self.pt.integrate(self.psit_nG, Q_ani, q=kpt.q)

        for Q_ni in Q_ani.values():
            self.assertTrue(Q_ni.dtype == self.dtype)

        P0_ani = dict([(a,Q_ni.copy()) for a,Q_ni in Q_ani.items()])
        self.pt.add(self.psit_nG, Q_ani, q=kpt.q)
        self.pt.integrate(self.psit_nG, P0_ani, q=kpt.q)

        #rank_a = self.gd.get_ranks_from_positions(spos_ac)
        #my_atom_indices = np.argwhere(self.gd.comm.rank == rank_a).ravel()

        #                                                ~ a   ~ a'
        #TODO XXX should fix PairOverlap-ish stuff for < p  | phi  > overlaps
        #                                                 i      i'

        #spos_ac = self.pt.spos_ac # NewLFC doesn't have this
        spos_ac = self.atoms.get_scaled_positions() % 1.0
        gpo = GridPairOverlap(self.gd, self.setups)
        B_aa = gpo.calculate_overlaps(spos_ac, self.pt)

        # Compare fingerprints across all processors
        fingerprint = np.array([md5_array(B_aa, numeric=True)])
        fingerprints = np.empty(world.size, np.int64)
        world.all_gather(fingerprint, fingerprints)
        if fingerprints.ptp(0).any():
            raise RuntimeError('Distributed matrices are not identical!')

        P_ani = dict([(a,Q_ni.copy()) for a,Q_ni in Q_ani.items()])
        for a1 in range(len(self.atoms)):
            if a1 in P_ani.keys():
                P_ni = P_ani[a1]
            else:
                # Atom a1 is not in domain so allocate a temporary buffer
                P_ni = np.zeros((self.bd.mynbands,self.setups[a1].ni,),
                                 dtype=self.dtype)
            for a2, Q_ni in Q_ani.items():
                # B_aa are the projector overlaps across atomic pairs
                B_ii = gpo.extract_atomic_pair_matrix(B_aa, a1, a2)
                P_ni += np.dot(Q_ni, B_ii.T) #sum over a2 and last i in B_ii
            self.gd.comm.sum(P_ni)

        self.check_and_plot(P_ani, P0_ani, 8, 'projection,linearity')

    def test_extrapolate_overlap(self):
        kpt = self.wfs.kpt_u[0]
        ppo = ProjectorPairOverlap(self.wfs, self.atoms)

        # Compare fingerprints across all processors
        fingerprint = np.array([md5_array(ppo.B_aa, numeric=True)])
        fingerprints = np.empty(world.size, np.int64)
        world.all_gather(fingerprint, fingerprints)
        if fingerprints.ptp(0).any():
            raise RuntimeError('Distributed matrices are not identical!')

        work_nG = np.empty_like(self.psit_nG)
        P_ani = ppo.apply(self.psit_nG, work_nG, self.wfs, kpt, \
            calculate_P_ani=True, extrapolate_P_ani=True)

        P0_ani = self.pt.dict(self.bd.mynbands)
        self.pt.integrate(work_nG, P0_ani, kpt.q)
        del work_nG

        self.check_and_plot(P_ani, P0_ani, 11, 'extrapolate,overlap')

    def test_extrapolate_inverse(self):
        kpt = self.wfs.kpt_u[0]
        ppo = ProjectorPairOverlap(self.wfs, self.atoms)

        # Compare fingerprints across all processors
        fingerprint = np.array([md5_array(ppo.B_aa, numeric=True)])
        fingerprints = np.empty(world.size, np.int64)
        world.all_gather(fingerprint, fingerprints)
        if fingerprints.ptp(0).any():
            raise RuntimeError('Distributed matrices are not identical!')

        work_nG = np.empty_like(self.psit_nG)
        P_ani = ppo.apply_inverse(self.psit_nG, work_nG, self.wfs, kpt, \
            calculate_P_ani=True, extrapolate_P_ani=True)

        P0_ani = self.pt.dict(self.bd.mynbands)
        self.pt.integrate(work_nG, P0_ani, kpt.q)
        del work_nG

        self.check_and_plot(P_ani, P0_ani, 11, 'extrapolate,inverse')

    def test_overlap_inverse_after(self):
        kpt = self.wfs.kpt_u[0]
        kpt.P_ani = self.pt.dict(self.bd.mynbands)
        ppo = ProjectorPairOverlap(self.wfs, self.atoms)

        # Compare fingerprints across all processors
        fingerprint = np.array([md5_array(ppo.B_aa, numeric=True)])
        fingerprints = np.empty(world.size, np.int64)
        world.all_gather(fingerprint, fingerprints)
        if fingerprints.ptp(0).any():
            raise RuntimeError('Distributed matrices are not identical!')

        work_nG = np.empty_like(self.psit_nG)
        self.pt.integrate(self.psit_nG, kpt.P_ani, kpt.q)
        P0_ani = dict([(a,P_ni.copy()) for a,P_ni in kpt.P_ani.items()])
        ppo.apply(self.psit_nG, work_nG, self.wfs, kpt, calculate_P_ani=False)

        res_nG = np.empty_like(self.psit_nG)
        ppo.apply_inverse(work_nG, res_nG, self.wfs, kpt, calculate_P_ani=True)
        del work_nG

        P_ani = self.pt.dict(self.bd.mynbands)
        self.pt.integrate(res_nG, P_ani, kpt.q)

        abserr = np.empty(1, dtype=float)
        for n in range(self.nbands):
            band_rank, myn = self.bd.who_has(n)
            if band_rank == self.bd.comm.rank:
                abserr[:] = np.abs(self.psit_nG[myn] - res_nG[myn]).max()
                self.gd.comm.max(abserr)
            self.bd.comm.broadcast(abserr, band_rank)
            self.assertAlmostEqual(abserr.item(), 0, 10)

        self.check_and_plot(P_ani, P0_ani, 10, 'overlap,inverse,after')

    def test_overlap_inverse_before(self):
        kpt = self.wfs.kpt_u[0]
        kpt.P_ani = self.pt.dict(self.bd.mynbands)
        ppo = ProjectorPairOverlap(self.wfs, self.atoms)

        # Compare fingerprints across all processors
        fingerprint = np.array([md5_array(ppo.B_aa, numeric=True)])
        fingerprints = np.empty(world.size, np.int64)
        world.all_gather(fingerprint, fingerprints)
        if fingerprints.ptp(0).any():
            raise RuntimeError('Distributed matrices are not identical!')

        work_nG = np.empty_like(self.psit_nG)
        self.pt.integrate(self.psit_nG, kpt.P_ani, kpt.q)
        P0_ani = dict([(a,P_ni.copy()) for a,P_ni in kpt.P_ani.items()])
        ppo.apply_inverse(self.psit_nG, work_nG, self.wfs, kpt, calculate_P_ani=False)

        res_nG = np.empty_like(self.psit_nG)
        ppo.apply(work_nG, res_nG, self.wfs, kpt, calculate_P_ani=True)
        del work_nG

        P_ani = self.pt.dict(self.bd.mynbands)
        self.pt.integrate(res_nG, P_ani, kpt.q)

        abserr = np.empty(1, dtype=float)
        for n in range(self.nbands):
            band_rank, myn = self.bd.who_has(n)
            if band_rank == self.bd.comm.rank:
                abserr[:] = np.abs(self.psit_nG[myn] - res_nG[myn]).max()
                self.gd.comm.max(abserr)
            self.bd.comm.broadcast(abserr, band_rank)
            self.assertAlmostEqual(abserr.item(), 0, 10)

        self.check_and_plot(P_ani, P0_ani, 10, 'overlap,inverse,before')

# -------------------------------------------------------------------

def UTGaussianWavefunctionFactory(boundaries, dtype):
    sep = '_'
    classname = 'UTGaussianWavefunctionSetup' \
    + sep + {'zero':'Zero', 'periodic':'Periodic', 'mixed':'Mixed'}[boundaries] \
    + sep + {float:'Float', complex:'Complex'}[dtype]
    class MetaPrototype(UTGaussianWavefunctionSetup, object):
        __doc__ = UTGaussianWavefunctionSetup.__doc__
        boundaries = boundaries
        dtype = dtype
    MetaPrototype.__name__ = classname
    return MetaPrototype

# -------------------------------------------------------------------

if __name__ in ['__main__', '__builtin__']:
    # We may have been imported by test.py, if so we should redirect to logfile
    if __name__ == '__builtin__':
        testrunner = CustomTextTestRunner('ut_invops.log', verbosity=2)
    else:
        stream = (world.rank == 0) and sys.stdout or devnull
        testrunner = TextTestRunner(stream=stream, verbosity=2)

    parinfo = []
    for test in [UTDomainParallelSetup_Zero, UTDomainParallelSetup_Periodic, \
                 UTDomainParallelSetup_Mixed]:
        info = ['', test.__name__, test.__doc__.strip('\n'), '']
        testsuite = initialTestLoader.loadTestsFromTestCase(test)
        map(testrunner.stream.writeln, info)
        testresult = testrunner.run(testsuite)
        assert testresult.wasSuccessful(), 'Initial verification failed!'
        parinfo.extend(['    Parallelization options: %s' % tci._parinfo for \
                        tci in testsuite._tests if hasattr(tci, '_parinfo')])
    parinfo = np.unique(np.sort(parinfo)).tolist()

    testcases = []
    for boundaries in ['zero', 'periodic', 'mixed']:
         for dtype in [float, complex]:
             testcases.append(UTGaussianWavefunctionFactory(boundaries, \
                 dtype))

    for test in testcases:
        info = ['', test.__name__, test.__doc__.strip('\n')] + parinfo + ['']
        testsuite = defaultTestLoader.loadTestsFromTestCase(test)
        map(testrunner.stream.writeln, info)
        testresult = testrunner.run(testsuite)
        # Provide feedback on failed tests if imported by test.py
        if __name__ == '__builtin__' and not testresult.wasSuccessful():
            raise SystemExit('Test failed. Check ut_invops.log for details.')

