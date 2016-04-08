
import sys
import time
import numpy as np

from gpaw.mpi import world
from gpaw.test.parunittest import ParallelTestCase, ParallelTextTestRunner, \
     defaultParallelTestLoader, main

# ------------------------------------------------------------------

class UTParallel(ParallelTestCase):
    """Parallel test suite of parallel test suite (don't ask)."""

    def reset_times(self, minwait=0.1, maxwait=0.2):
        if world.size == 1:
            return
        self.time_r = np.random.uniform(minwait, maxwait, size=world.size)
        world.broadcast(self.time_r, 0)

    def sleep_times(self):
        if world.size == 1:
            return
        time.sleep(self.time_r[world.rank])


class UTParallel_Succeeds(UTParallel):
    __doc__ = UTParallel.__doc__
    passes = True
    numtests = 1
    numerrors = 0
    numfails = 0

    def test_all_passes(self):
        self.reset_times()
        self.assertTrue(True)
        self.sleep_times()


class UTParallel_Raises(UTParallel):
    __doc__ = UTParallel.__doc__
    passes = False
    numtests = 2
    numerrors = 2
    numfails = 0

    def test_master_raises(self):
        self.reset_times()
        if world.rank == 0:
            raise Exception('I am special!')
        self.sleep_times()

    # Nichols A. Romero
    # naromero@alcf.anl.gov
    # Feb. 28, 2010
    # This test may not be valid
    # def test_slave_exits(self):
    #     self.reset_times()
    #     if world.rank == world.size-1:
    #         sys.exit(-1)
    #     self.sleep_times()

    def test_some_raises(self):
        self.reset_times()
        if world.rank % 3 == 0:
            raise RuntimeError
        elif world.rank % 3 == 1:
            raise IOError
        self.sleep_times()


class UTParallel_Fails(UTParallel):
    __doc__ = UTParallel.__doc__
    passes = False
    numtests = 4
    numerrors = 0
    numfails = 4

    def test_master_fails(self):
        self.reset_times()
        self.assertTrue(world.rank!=0)
        self.sleep_times()

    def test_all_fails(self):
        self.reset_times()
        self.assertTrue(False)
        self.sleep_times()

    def test_odd_fails(self):
        self.reset_times()
        self.assertTrue(world.size>1 and world.rank%2==0)
        self.sleep_times()

    def test_even_fails(self):
        self.reset_times()
        self.assertTrue(world.size>1 and world.rank%2==1)
        self.sleep_times()


# ------------------------------------------------------------------

if __name__ in ['__main__', '__builtin__']:
    # We may have been imported by test.py, if so we should redirect to logfile
    output = (__name__ == '__builtin__') and 'ut_parallel.log' or sys.stdout
    testrunner = ParallelTextTestRunner(stream=output, verbosity=2)

    testcases = [UTParallel_Succeeds, UTParallel_Raises, UTParallel_Fails]

    for test in testcases:
        info = '\n' + test.__name__ + '\n' + test.__doc__.strip('\n') + '\n'
        testsuite = defaultParallelTestLoader.loadTestsFromTestCase(test)
        testrunner.stream.writeln(info)
        testresult = testrunner.run(testsuite)
        verified = ((testresult.wasSuccessful() == test.passes) and
            (testresult.testsRun == test.numtests) and
            (len(testresult.errors) == test.numerrors) and
            (len(testresult.failures) == test.numfails))
        if verified:
            conclusion = 'STATUS: %s performed as expected.' % test
        else:
            conclusion = 'STATUS: %s did not meet expectations.\n\n' \
                ' ** Number of tests:    %2d ran,    %2d expected.\n%s\n' \
                ' ** Number of errors:   %2d raised, %2d expected.\n%s\n' \
                ' ** Number of failures: %2d failed, %2d expected.\n%s' \
                % (test, testresult.testsRun, test.numtests, testsuite, \
                   len(testresult.errors), test.numerrors, testresult.errors, \
                   len(testresult.failures), test.numfails, testresult.failures)

        separator0 = '*'*len(testresult.separator1)
        map(testrunner.stream.writeln, ['',separator0,conclusion,separator0,''])

        # Provide feedback on failed tests if imported by test.py
        if __name__ == '__builtin__' and not verified:
            raise SystemExit('Test failed. Check ut_parallel.log for details.')

