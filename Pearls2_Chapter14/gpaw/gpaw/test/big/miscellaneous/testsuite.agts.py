def agts(queue):
    queue.add('testsuite.agts.py', ncpus=8, walltime=80)

if __name__ == '__main__':
    # Run test suite
    from gpaw.test.test import run
    nfailed = run()
    assert nfailed == 0
