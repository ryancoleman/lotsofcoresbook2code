.. _testing:

============
Testing GPAW
============

Testing of gpaw is done by a nightly test suite consisting of many
small and quick tests and by a weekly set of larger test.


Quick test suite
================

Use the :program:`gpaw-test` command to run the tests::

    $ gpaw-test --help
    Usage: gpaw-test [options] [tests]
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -x test1.py,test2.py,..., --exclude=test1.py,test2.py,...
                            Exclude tests (comma separated list of tests).
      -f, --run-failed-tests-only
                            Run failed tests only.
      --from=TESTFILE       Run remaining tests, starting from TESTFILE
      --after=TESTFILE      Run remaining tests, starting after TESTFILE
      --range=test_i.py,test_j.py
                            Run tests in range test_i.py to test_j.py (inclusive)
      -j JOBS, --jobs=JOBS  Run JOBS threads.
      --reverse             Run tests in reverse order (less overhead with
                            multiple jobs)
      -k, --keep-temp-dir   Do not delete temporary files.

A temporary directory will be made and the tests will run in that
directory.  If all tests pass, the directory is removed.

The test suite consists of a large number of small and quick tests
found in the :trac:`gpaw/test` directory.  The tests run nightly in serial
and in parallel.  Here are the results from `BuildBot
<https://gpaw-buildbot.fysik.dtu.dk/waterfall>`_.



Adding new tests
----------------

A test script should fulfill a number of requirements:

* It should be quick.  Preferably a few seconds, but a few minutes is
  OK.  If the test takes several minutes or more, consider making the
  test a :ref:`big test <big-test>`.

* It should not depend on other scripts.

* It should be possible to run it on 1, 2, 4, and 8 cores.

A test can produce standard output and files - it doesn't have to
clean up.  Remember to add the new test to list of all tests specified
in the :trac:`gpaw/test/__init__.py` file.

Use this function to check results:

.. function:: gpaw.test.equal(x, y, tolerance=0, fail=True, msg='')


.. _big-test:

Big tests
=========

The directory in :trac:`gpaw/test/big` contains a set of longer and more
realistic tests that we run every weekend.  These are submitted to a
queueing system of a large computer.  Here is an example for Niflheim:
:trac:`gpaw/test/big/niflheim.py`.


Adding new tests
----------------

To add a new test, create a script somewhere in the file hierarchy ending with
.agts.py (e.g. ``submit.agts.py``). ``AGTS`` is short for Advanced GPAW Test
System (or Another Great Time Sink). This script defines how a number of
scripts should be submitted to niflheim and how they depend on each other.
Consider an example where one script, calculate.py, calculates something and
saves a .gpw file and another script, analyse.py, analyses this output. Then
the submit script should look something like::

    def agts(queue):
        calc = queue.add('calculate.py',
                         ncpus=8,
                         walltime=25)

        queue.add('analyse.py'
                  ncpus=1,
                  walltime=5,
                  deps=[calc])

As shown, this script has to contain the definition of the function
``agts()`` which should take exactly one argument, ``queue``. Then
each script is added to a queue object, along with some data which
defines how and when to run the job.  Note how ``queue.add()`` returns
a job object which can be used to specify dependencies.

Possible keys are:

=============  ========  =============  ===================================
Name           Type      Default value  Description
=============  ========  =============  ===================================
``ncpus``      ``int``   ``1``          Number of cpus
``walltime``   ``int``   ``15``         Requested walltime in minutes
``deps``       ``list``  ``[]``         List of jobs this job depends on
``creates``    ``list``  ``[]``         List of files this job creates
                                        (figures and other stuff for the
                                        web-page) 
=============  ========  =============  ===================================


Creating expensive figures for the web-page
-------------------------------------------

You can also use AGTS for running long jobs that create figures for this
web-page.  For an example, look at the source code
:trac:`here <doc/tutorials/stm>` which will produce this: :ref:`stm tutorial`.
