from __future__ import print_function
import os, sys, time
import numpy as np

try:
    # Matplotlib is not a dependency
    import matplotlib as mpl
    mpl.use('Agg')  # force the antigrain backend
except (ImportError, RuntimeError):
    mpl = None

from ase import Atoms
from ase.structure import molecule
from ase.parallel import paropen
from ase.units import Bohr, Hartree
from ase.io.trajectory import PickleTrajectory
from ase.calculators.singlepoint import SinglePointCalculator
from ase.utils import devnull

from gpaw import GPAW, debug
from gpaw.mpi import world
from gpaw.tddft import TDDFT
from gpaw.tddft.units import attosec_to_autime

# -------------------------------------------------------------------

from gpaw.test.ut_common import ase_svnversion, shapeopt, TestCase, \
    TextTestRunner, CustomTextTestRunner, defaultTestLoader, \
    initialTestLoader, create_random_atoms, create_parsize_maxbands

# -------------------------------------------------------------------

class UTGroundStateSetup(TestCase):
    """
    Setup a simple calculation starting from ground state."""

    # Number of bands
    nbands = 1

    # Spin-paired, single kpoint
    nspins = 1
    nibzkpts = 1

    # Mean spacing and number of grid points per axis (G x G x G)
    h = 0.25 / Bohr

    def setUp(self):
        for virtvar in []:
            assert getattr(self,virtvar) is not None, 'Virtual "%s"!' % virtvar

        parsize_domain, parsize_bands = create_parsize_maxbands(self.nbands, world.size)
        self.parallel = {'domain': parsize_domain, 'band': parsize_bands}

        self.atoms = molecule('Na2')
        self.atoms.center(vacuum=4.0)
        self.atoms.set_pbc(False)

        cell_c = np.sum(self.atoms.get_cell()**2, axis=1)**0.5 / Bohr
        ngpts = 16 * np.round(cell_c / (self.h * 16))
        self.gsname = 'ut_tddft_gs'
        self.gscalc = GPAW(gpts=ngpts, nbands=self.nbands, basis='dzp',
                           setups={'Na': '1'},
                           spinpol=(self.nspins==2), parallel=self.parallel,
                           txt=self.gsname + '.txt')

    def tearDown(self):
        del self.atoms, self.gscalc

    # =================================

    def verify_comm_sizes(self):
        if world.size == 1:
            return
        rem = world.size // (self.parallel['band'] * self.parallel['domain'])
        comm_sizes = (world.size, self.parallel['band'],
                      self.parallel['domain'], rem)
        self._parinfo =  '%d world, %d band, %d domain, %d kpt' % comm_sizes
        self.assertEqual(self.nbands % self.parallel['band'], 0)
        self.assertEqual((self.nspins * self.nibzkpts) % rem, 0)

    def verify_ground_state(self):
        #XXX DEBUG START
        if debug and os.path.isfile(self.gsname + '.gpw'):
            return
        #XXX DEBUG END

        self.atoms.set_calculator(self.gscalc)
        self.assertAlmostEqual(self.atoms.get_potential_energy(), -1.0621, 4)
        self.gscalc.write(self.gsname + '.gpw', mode='all')
        self.assertTrue(os.path.isfile(self.gsname + '.gpw'))

# -------------------------------------------------------------------

class UTStaticPropagatorSetup(UTGroundStateSetup):
    __doc__ = UTGroundStateSetup.__doc__ + """
    Propagate electrons with fixed nuclei and verify results.""" #TODO

    duration = 10.0 #24.0

    timesteps = None
    propagator = None

    def setUp(self):
        UTGroundStateSetup.setUp(self)

        for virtvar in ['timesteps', 'propagator']:
            assert getattr(self,virtvar) is not None, 'Virtual "%s"!' % virtvar

        self.tdname = 'ut_tddft_td_' + self.propagator.lower()
        self.tdcalc = TDDFT(self.gsname + '.gpw', propagator=self.propagator,
                            txt=self.tdname + '.txt')

    def tearDown(self):
        del self.tdcalc

    # =================================

    def _test_timestepping(self, t):
        #XXX DEBUG START
        if debug and os.path.isfile('%s_%d.gpw' % (self.tdname, t)):
            return
        #XXX DEBUG END

        timestep = self.timesteps[t]
        self.assertAlmostEqual(self.duration % timestep, 0.0, 12)
        niter = int(self.duration / timestep)
        ndiv = 1 #XXX

        traj = PickleTrajectory('%s_%d.traj' % (self.tdname, t),
                                'w', self.tdcalc.get_atoms())

        t0 = time.time()
        f = paropen('%s_%d.log' % (self.tdname, t), 'w')
        print('propagator: %s, duration: %6.1f as, timestep: %5.2f as, ' \
            'niter: %d' % (self.propagator, self.duration, timestep, niter), file=f)

        for i in range(1, niter+1):
            # XXX bare bones propagation without all the nonsense
            self.tdcalc.propagator.propagate(self.tdcalc.time,
                                             timestep * attosec_to_autime)
            self.tdcalc.time += timestep * attosec_to_autime
            self.tdcalc.niter += 1

            if i % ndiv == 0:
                rate = 60 * ndiv / (time.time()-t0)
                ekin = self.tdcalc.atoms.get_kinetic_energy()
                epot = self.tdcalc.get_td_energy() * Hartree
                F_av = np.zeros((len(self.tdcalc.atoms), 3))
                print('i=%06d, time=%6.1f as, rate=%6.2f min^-1, ' \
                    'ekin=%13.9f eV, epot=%13.9f eV, etot=%13.9f eV' \
                    % (i, timestep * i, rate, ekin, epot, ekin + epot), file=f)
                t0 = time.time()

                # Hack to prevent calls to GPAW::get_potential_energy when saving
                spa = self.tdcalc.get_atoms()
                spc = SinglePointCalculator(epot, F_av, None, None, spa)
                spa.set_calculator(spc)
                traj.write(spa)
        f.close()
        traj.close()
        self.tdcalc.write('%s_%d.gpw' % (self.tdname, t), mode='all')

        # Save density and wavefunctions to binary
        gd, finegd = self.tdcalc.wfs.gd, self.tdcalc.density.finegd
        if world.rank == 0:
            big_nt_g = finegd.collect(self.tdcalc.density.nt_g)
            np.save('%s_%d_nt.npy' % (self.tdname, t), big_nt_g)
            del big_nt_g

            big_psit_nG = gd.collect(self.tdcalc.wfs.kpt_u[0].psit_nG)
            np.save('%s_%d_psit.npy' % (self.tdname, t), big_psit_nG)
            del big_psit_nG
        else:
            finegd.collect(self.tdcalc.density.nt_g)
            gd.collect(self.tdcalc.wfs.kpt_u[0].psit_nG)
        world.barrier()

    def test_timestepping_ref(self):
        # Reference values for density and wavefunctions (assumed stationary)
        nt0_g = self.tdcalc.density.nt_g
        phases_n = self.tdcalc.wfs.kpt_u[0].eps_n * self.duration * attosec_to_autime
        psit0_nG = np.exp(-1j * phases_n) * self.tdcalc.wfs.kpt_u[0].psit_nG

        f = paropen('%s_ref.log' % self.tdname, 'w')
        niters = np.round(self.duration / self.timesteps).astype(int)
        print('propagator: %s, duration: %6.1f as, niters: %s, ' \
            % (self.propagator, self.duration, niters.tolist()), file=f)

        for t,timestep in enumerate(self.timesteps):
            self.assertTrue(os.path.isfile('%s_%d.gpw' % (self.tdname, t)))

            # Load density and wavefunctions from binary
            gd, finegd = self.tdcalc.wfs.gd, self.tdcalc.density.finegd
            nt_g, psit_nG = finegd.empty(), gd.empty(self.nbands, dtype=complex)
            if world.rank == 0:
                big_nt_g = np.load('%s_%d_nt.npy' % (self.tdname, t))
                finegd.distribute(big_nt_g, nt_g)
                del big_nt_g

                big_psit_nG = np.load('%s_%d_psit.npy' % (self.tdname, t))
                gd.distribute(big_psit_nG, psit_nG)
                del big_psit_nG
            else:
                finegd.distribute(None, nt_g)
                gd.distribute(None, psit_nG)

            # Check loaded density and wavefunctions against reference values
            dnt = finegd.comm.max(np.abs(nt_g - nt0_g).max())
            dpsit = gd.comm.max(np.abs(psit_nG - psit0_nG).max())
            print('t=%d, timestep: %5.2f as, dnt: %16.13f, ' \
                'dpsit: %16.13f' % (t, timestep, dnt, dpsit), file=f)
            snt, spsit = {'SITE': (5,4)}.get(self.propagator, (7,5))
            #self.assertAlmostEqual(dnt, 0, snt, 't=%d, timestep: ' \
            #    '%5.2f as, dnt: %g, digits: %d' % (t, timestep, dnt, snt))
            #self.assertAlmostEqual(dpsit, 0, spsit, 't=%d, timestep: ' \
            #    '%5.2f as, dpsit=%g, digits: %d' % (t, timestep, dpsit, spsit))
        f.close()

# -------------------------------------------------------------------

import new

def UTStaticPropagatorFactory(timesteps, propagator):
    sep = '_'
    classname = 'UTStaticPropagatorSetup' \
    + sep + propagator
    class MetaPrototype(UTStaticPropagatorSetup, object):
        __doc__ = UTStaticPropagatorSetup.__doc__
        timesteps = timesteps
        propagator = propagator
    for t,timestep in enumerate(timesteps):
        func = lambda _self, _t=t: _self._test_timestepping(_t)
        method = new.instancemethod(func, None, MetaPrototype)
        method_name = 'test_timestepping_%02.0fas' % timestep
        setattr(MetaPrototype, method_name, method)
    MetaPrototype.__name__ = classname
    return MetaPrototype

# -------------------------------------------------------------------

if __name__ in ['__main__', '__builtin__']:
    # We may have been imported by test.py, if so we should redirect to logfile
    if __name__ == '__builtin__':
        testrunner = CustomTextTestRunner('ut_tddft.log', verbosity=2)
    else:
        stream = (world.rank == 0) and sys.stdout or devnull
        testrunner = TextTestRunner(stream=stream, verbosity=2)

    parinfo = []
    for test in [UTGroundStateSetup]:
        info = ['', test.__name__, test.__doc__.strip('\n'), '']
        testsuite = initialTestLoader.loadTestsFromTestCase(test)
        map(testrunner.stream.writeln, info)
        testresult = testrunner.run(testsuite)
        assert testresult.wasSuccessful(), 'Initial verification failed!'
        parinfo.extend(['    Parallelization options: %s' % tci._parinfo for \
                        tci in testsuite._tests if hasattr(tci, '_parinfo')])
    parinfo = np.unique(np.sort(parinfo)).tolist()

    # Try all timesteps between 1 and 10 that factor the duration nicely
    timesteps = np.arange(1, min(10, UTStaticPropagatorSetup.duration) + 1)
    timesteps = timesteps[UTStaticPropagatorSetup.duration % timesteps == 0]

    testcases = []
    for propagator in ['ECN', 'SICN', 'SITE', 'SIKE', 'ETRSCN']:
        testcases.append(UTStaticPropagatorFactory(timesteps, propagator))

    for test in testcases:
        info = ['', test.__name__, test.__doc__.strip('\n')] + parinfo + ['']
        testsuite = defaultTestLoader.loadTestsFromTestCase(test)
        map(testrunner.stream.writeln, info)
        testresult = testrunner.run(testsuite)
        # Provide feedback on failed tests if imported by test.py
        if __name__ == '__builtin__' and not testresult.wasSuccessful():
            raise SystemExit('Test failed. Check ut_tddft.log for details.')

