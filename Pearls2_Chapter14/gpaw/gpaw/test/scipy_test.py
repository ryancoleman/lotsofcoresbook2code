from __future__ import print_function
import sys
from scipy import test

from gpaw.mpi import rank

_stdout = sys.stdout
_stderr = sys.stderr
# scipy tests write to stderr
sys.stderr = open("scipy_test%02d.out" % rank, "w")
result = test(verbose=10)
sys.stdout = _stdout
sys.stderr = _stderr
if not result.wasSuccessful():
    print("scipy_test%02d.out" % rank, result.errors, result.failures, file=sys.stderr)
assert result.wasSuccessful()

