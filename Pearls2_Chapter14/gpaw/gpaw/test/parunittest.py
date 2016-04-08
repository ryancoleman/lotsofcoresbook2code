"""
Python unit testing framework, based on Erich Gamma's JUnit and Kent Beck's
Smalltalk testing framework. Parallel support by Christian Glinsvad.

For further information please consult the documentation for the Unittest
module, from which this module is derived. Online reference manual at:

  http://docs.python.org/dev/library/unittest.html

Copyright (c) 1999-2003 Steve Purcell
This module is free software, and you may redistribute it and/or modify
it under the same terms as Python itself, so long as this copyright message
and disclaimer are retained in their original form.

IN NO EVENT SHALL THE AUTHOR BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF
THIS CODE, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

THE AUTHOR SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE.  THE CODE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS,
AND THERE IS NO OBLIGATION WHATSOEVER TO PROVIDE MAINTENANCE,
SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
"""
from __future__ import print_function
import sys
import numpy as np
from unittest import TestResult, TestCase, TestSuite, \
                     _TextTestResult, TextTestRunner, TestLoader, \
                     FunctionTestCase, TestProgram, defaultTestLoader

from ase.utils import devnull

from gpaw import debug
from gpaw.mpi import world, broadcast_string

# -------------------------------------------------------------------
# Exported classes and functions
# -------------------------------------------------------------------

__all__ = ['ParallelTestResult', 'ParallelTestCase', 'ParallelTestSuite',
           'ParallelTextTestRunner', 'ParallelTestLoader',
           'ParallelFunctionTestCase', 'main', 'defaultParallelTestLoader']

# -------------------------------------------------------------------
# Limited backward compatibility
# -------------------------------------------------------------------

if sys.version_info[:2] < (2, 3):
    raise RuntimeError('Python 2.3 or greater required!')

try:
    from unittest import __version__
except ImportError:
    unittest_version = (2, 0)
else:
    unittest_version = tuple(map(int, __version__.split('.')))

# User interface should at least comply with Unittest version 1.56 rev. 34209
if unittest_version < (1,56):
    _TestCase = TestCase
    class TestCase(_TestCase):
        assertTrue = _TestCase.failUnless
        assertFalse = _TestCase.failIf

# -------------------------------------------------------------------
# Test framework core
# -------------------------------------------------------------------

class ParallelTestResult(TestResult):
    __doc__ = TestResult.__doc__

    def __init__(self, comm=None):
        if comm is None:
            comm = world
        self.comm = comm
        self.outcomes = []
        self.last_errors = np.empty(self.comm.size, dtype=bool)
        self.last_failed = np.empty(self.comm.size, dtype=bool)
        TestResult.__init__(self)

    def startTest(self, test):
        "Called when the given test is about to be run. Resets global flags."
        self.last_errors.fill(False)
        self.last_failed.fill(False)
        TestResult.startTest(self, test)

    def stopTest(self, test):
        """Called when the given test has been run. If the stop flag was
        raised beforehand, will broadcast to raise flags for global stop."""
        stop_flags = np.empty(self.comm.size, dtype=bool)
        self.comm.all_gather(np.array([self.shouldStop]), stop_flags)
        self.shouldStop = stop_flags.any()
        TestResult.stopTest(self, test)

    def _exc_info_to_string(self, err, test):
        # Includes second argument as of Unittest version 1.63 rev. 34859
        if unittest_version < (1,63):
            return TestResult._exc_info_to_string(self, err)
        else:
            return TestResult._exc_info_to_string(self, err, test)

    def addStatus(self, test, err=None, error=False, failure=False):
        """Negotiate global status immediately after a test has been run.
        Called on all processor with the local test status (i.e. error, failure
        or success). 'err' is a tuple of values as returned by sys.exc_info().
        """
        if error and failure:
            raise RuntimeError('Parallel unittest can\'t handle simultaneous' \
                               + ' errors and failures within a single test.')

        self.comm.all_gather(np.array([error]), self.last_errors)
        if self.last_errors.any():
            all_texts = []
            for rank in np.argwhere(self.last_errors).ravel():
                if rank == self.comm.rank:
                    assert self.last_errors[self.comm.rank]
                    text = self._exc_info_to_string(err, test)
                else:
                    text = None
                text = broadcast_string(text, root=rank, comm=self.comm)
                all_texts.append((rank,text))
            self.errors.append((test, all_texts))

        self.comm.all_gather(np.array([failure]), self.last_failed)
        if self.last_failed.any():
            all_texts = []
            for rank in np.argwhere(self.last_failed).ravel():
                if rank == self.comm.rank:
                    assert self.last_failed[self.comm.rank]
                    text = self._exc_info_to_string(err, test)
                else:
                    text = None
                text = broadcast_string(text, root=rank, comm=self.comm)
                all_texts.append((rank,text))
            self.failures.append((test, all_texts))

    def addError(self, test, err):
        """Negotiate global status. Called when an error has occurred on a
        processor. 'err' is a tuple of values as returned by sys.exc_info().
        """
        self.addStatus(test, err, error=True)

    def addFailure(self, test, err):
        """Negotiate global status. Called when an failure has occurred on a
        processor. 'err' is a tuple of values as returned by sys.exc_info().
        """
        self.addStatus(test, err, failure=True)

    def addSuccess(self, test):
        """Probe global status for potential errors by passive negotiation.
        Called when a test has completed successfully on a processor.
        """
        self.addStatus(test)

    def addSkip(self, test, reason):
        """Called when a test is skipped. Not ready!"""
        raise NotImplementedError

    def addExpectedFailure(self, test, err):
        """Called when an expected failure/error occured. Not ready!"""
        raise NotImplementedError

    def addUnexpectedSuccess(self, test):
        """Called when a test was expected to fail, but succeed. Not ready!"""
        raise NotImplementedError


class ParallelTestCase(TestCase):
    __doc__ = TestCase.__doc__

    def defaultTestResult(self):
        return ParallelTestResult()

    def debug(self):
        """Run the test without collecting errors in a TestResult"""
        print('WARNING: Test case is strictly serial in debug mode!')
        TestCase.debug(self)


# No point in implementing these as the originals work just fine
ParallelTestSuite = TestSuite
ParallelFunctionTestCase = FunctionTestCase


# -------------------------------------------------------------------
# Locating and loading tests
# -------------------------------------------------------------------

# No point in implementing these as the originals work just fine
ParallelTestLoader = TestLoader
defaultParallelTestLoader = ParallelTestLoader()


# -------------------------------------------------------------------
# Text UI
# -------------------------------------------------------------------

class _ParallelTextTestResult(ParallelTestResult, _TextTestResult):
    __doc__ = _TextTestResult.__doc__
    allErrors = False

    def __init__(self, comm, *args, **kwargs):
        ParallelTestResult.__init__(self, comm)
        _TextTestResult.__init__(self, *args, **kwargs) # use new-style?

    def startTest(self, test):
        ParallelTestResult.startTest(self, test)
        if self.showAll:
            self.stream.write(self.getDescription(test))
            self.stream.write(" ... ")
            self.stream.flush()

    def addSuccess(self, test):
        ParallelTestResult.addSuccess(self, test)

    def addError(self, test, err):
        ParallelTestResult.addError(self, test, err)

    def addFailure(self, test, err):
        ParallelTestResult.addFailure(self, test, err)

    def shortRanks(self, ranks):
        if len(ranks) == 0:
            return 'none'
        elif isinstance(ranks, np.ndarray) and ranks.dtype == bool:
            if ranks.all():
                return 'all'
            ranks = np.argwhere(ranks).ravel()
        ranks = np.sort(ranks)
        if np.all(ranks == range(self.comm.size)):
            return 'all'
        return 'rank ' + ','.join(map(str,ranks))

    def stopTest(self, test):
        self.stream.flush()
        self.comm.barrier()

        if self.last_errors.any() and self.last_failed.any():
            if not debug:
                raise RuntimeError('Parallel unittest can\'t handle ' \
                    'simultaneous errors and failures within a single test.')
            if self.showAll:
                error_rankinfo = self.shortRanks(self.last_errors)
                fail_rankinfo = self.shortRanks(self.last_failed)
                self.stream.writeln("ERROR (%s) / FAIL (%s)" \
                                    % (error_rankinfo, fail_rankinfo))
            elif self.dots:
                self.stream.writeln('M')
        elif self.last_errors.any():
            if self.showAll:
                rankinfo = self.shortRanks(self.last_errors)
                self.stream.writeln("ERROR (%s)" % rankinfo)
            elif self.dots:
                self.stream.writeln('E')
        elif self.last_failed.any():
            if self.showAll:
                rankinfo = self.shortRanks(self.last_failed)
                self.stream.writeln("FAIL (%s)" % rankinfo)
            elif self.dots:
                self.stream.writeln('F')
        else:
            if self.showAll:
                self.stream.writeln("ok")
            else:
                self.stream.writeln('.')

    def printErrorList(self, flavour, errors):
        for test, all_texts in errors:
            description = self.getDescription(test)
            if all_texts:
                self.stream.writeln(self.separator1)
            if self.allErrors:
                for rank, text2 in all_texts:
                    rankinfo = self.shortRanks([rank]).upper()
                    text1 = '%s-%s: %s' % (flavour, rankinfo, description)
                    self.stream.writeln(text1)
                    self.stream.writeln(self.separator2)
                    self.stream.writeln(text2)
            else:
                text_ranks = {}
                for rank, text2 in all_texts:
                    if text2 not in text_ranks:
                        text_ranks[text2] = []
                    text_ranks[text2].append(rank)
                for text2, ranks in text_ranks.items():
                    rankinfo = self.shortRanks(ranks).upper()
                    text1 = '%s-%s: %s' % (flavour, rankinfo, description)
                    self.stream.writeln(text1)
                    self.stream.writeln(self.separator2)
                    self.stream.writeln(text2)


class ParallelTextTestRunner(TextTestRunner):
    __doc__ = TextTestRunner.__doc__
    logfile = None

    def __init__(self, comm=None, stream=sys.stderr, **kwargs):
        if comm is None:
            comm = world
        self.comm = comm
        if self.comm.rank != 0:
            stream = devnull
        elif type(stream) is str:
            self.logfile = stream
            stream = open(self.logfile, 'w', buffering=0)
        TextTestRunner.__init__(self, stream=stream, **kwargs)

    def _makeResult(self):
        return _ParallelTextTestResult(self.comm, self.stream, \
            self.descriptions, self.verbosity)

    def run(self, test):
        stderr_old = sys.stderr
        try:
            sys.stderr = self.stream
            testresult = TextTestRunner.run(self, test)
        finally:
            sys.stderr = stderr_old
        return testresult


# -------------------------------------------------------------------
# Facilities for running tests from the command line
# -------------------------------------------------------------------

class ParallelTestProgram(TestProgram):
    __doc__ = TestProgram.__doc__

    def runTests(self):
        if self.testRunner is None:
            self.testRunner = ParallelTextTestRunner(verbosity=self.verbosity)
        TestProgram.runTests(self)

main = ParallelTestProgram

