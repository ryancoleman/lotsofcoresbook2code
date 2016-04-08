# Submit tests from the test suite that were remove becase they were
# too long.
def agts(queue):
    queue.add('H2Al110.py')
    queue.add('dscf_CO.py')
    queue.add('revtpss_tpss_scf.py')
    queue.add('ltt.py')
    queue.add('pblacs_oblong.py', walltime=5, ncpus=64)
