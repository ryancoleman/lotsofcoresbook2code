from __future__ import print_function
import os
import platform
import tempfile
import warnings
from optparse import OptionParser

import gpaw.mpi as mpi
from gpaw import debug
from gpaw.version import version


def main(args=None):
    description = ('Run the GPAW test suite.  The test suite can be run in '
                   'parallel with MPI through gpaw-python.  The test suite '
                   'supports 1, 2, 4 or 8 CPUs although some tests are '
                   'skipped for some parallelizations.  If no TESTs are '
                   'given, run all tests supporting the parallelization.')

    parser = OptionParser(usage='%prog [OPTION...] [TEST...]',
                          description=description,
                          version='%%prog %s' % version)
    parser.add_option('-x', '--exclude',
                      type='string', default=None,
                      help='Exclude tests (comma separated list of tests).',
                      metavar='test1.py,test2.py,...')
    parser.add_option('-f', '--run-failed-tests-only',
                      action='store_true',
                      help='Run failed tests only.')
    parser.add_option('--from', metavar='TESTFILE', dest='from_test',
                      help='Run remaining tests, starting from TESTFILE')
    parser.add_option('--after', metavar='TESTFILE', dest='after_test',
                      help='Run remaining tests, starting after TESTFILE')
    parser.add_option('--range',
                      type='string', default=None,
                      help='Run tests in range test_i.py to test_j.py '
                      '(inclusive)',
                      metavar='test_i.py,test_j.py')
    parser.add_option('-j', '--jobs', type='int', default=1,
                      help='Run JOBS threads.  Each test will be executed '
                      'in serial by one thread.  This option cannot be used '
                      'for parallelization together with MPI.')
    parser.add_option('--reverse', action='store_true',
                      help=('Run tests in reverse order (less overhead with '
                            'multiple jobs)'))
    parser.add_option('-k', '--keep-temp-dir', action='store_true',
                      dest='keep_tmpdir',
                      help='Do not delete temporary files.')
    parser.add_option('-d', '--directory', help='Run test in this directory')
    parser.add_option('-s', '--show-output', action='store_true',
                      help='Show standard output from tests.')

    opt, tests = parser.parse_args(args)

    if len(tests) == 0:
        from gpaw.test import tests

    if opt.reverse:
        tests.reverse()

    if opt.run_failed_tests_only:
        tests = [line.strip() for line in open('failed-tests.txt')]

    exclude = []
    if opt.exclude is not None:
        exclude += opt.exclude.split(',')

    if opt.from_test:
        fromindex = tests.index(opt.from_test)
        tests = tests[fromindex:]

    if opt.after_test:
        index = tests.index(opt.after_test) + 1
        tests = tests[index:]

    if opt.range:
        # default start(stop) index is first(last) test
        indices = opt.range.split(',')
        try:
            start_index = tests.index(indices[0])
        except ValueError:
            start_index = 0
        try:
            stop_index = tests.index(indices[1]) + 1
        except ValueError:
            stop_index = len(tests)
        tests = tests[start_index:stop_index]

    if opt.jobs > 1:
        exclude.append('maxrss.py')

    for test in exclude:
        if test in tests:
            tests.remove(test)

    from gpaw.test import TestRunner

    if mpi.world.size > 8:
        if mpi.rank == 0:
            message = '!!!!!!!\n' \
                'GPAW regression test suite was not designed to run on more\n' \
                'than 8 MPI tasks. Re-run test suite using 1, 2, 4 or 8 MPI\n' \
                'tasks instead.'
            warnings.warn(message, RuntimeWarning)

    if mpi.rank == 0:
        if opt.directory is None:
            tmpdir = tempfile.mkdtemp(prefix='gpaw-test-')
        else:
            tmpdir = opt.directory
            if os.path.isdir(tmpdir):
                opt.keep_tmpdir = True
            else:
                os.mkdir(tmpdir)
    else:
        tmpdir = None
    tmpdir = mpi.broadcast_string(tmpdir)
    cwd = os.getcwd()
    os.chdir(tmpdir)
    operating_system = platform.system() + ' ' + platform.machine()
    operating_system += ' ' + ' '.join(platform.dist())
    python = platform.python_version() + ' ' + platform.python_compiler()
    python += ' ' + ' '.join(platform.architecture())
    if mpi.rank == 0:
        print('python %s on %s' % (python, operating_system))
        print('Running tests in %s' % tmpdir)
        print('Jobs: %d, Cores: %d, debug-mode: %r' % (opt.jobs, mpi.size,
                                                       debug))
    failed = TestRunner(tests, jobs=opt.jobs,
                        show_output=opt.show_output).run()
    os.chdir(cwd)
    if mpi.rank == 0:
        if len(failed) > 0:
            open('failed-tests.txt', 'w').write('\n'.join(failed) + '\n')
        elif not opt.keep_tmpdir:
            os.system('rm -rf ' + tmpdir)
    return len(failed)

# Backwards compatibility:
run = main

if __name__ == '__main__':
    main()
