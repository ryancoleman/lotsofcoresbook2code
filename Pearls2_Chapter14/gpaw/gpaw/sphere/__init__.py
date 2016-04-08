
import numpy as np

from gpaw.utilities import ffact

# Define (l+|m|)!/(l-|m|)!
lmfact = lambda l,m: ffact(l-abs(m), l+abs(m))

def _lmiter(lmax, full=True):
    for l in xrange(lmax+1):
        for m in xrange(full and -l or 0,l+1):
            yield (l,m)
    raise StopIteration

def lmiter(lmax, full=True, comm=None, cost=None):
    """Utility function to parallelize over (l,m) with load balancing."""
    if comm is None or comm.size == 1:
        return _lmiter(lmax ,full)

    lm_j = tuple(_lmiter(lmax, full))
    if cost is None:
        cost_j = np.ones(len(lm_j), dtype=float)
    else:
        cost_j = np.fromiter((cost(lmax, *lm) for lm in lm_j), dtype=float)
    assert np.isfinite(cost_j).all() and np.all(cost_j>0) #XXX and np.isreal?
    a,b =  np.array([comm.rank,comm.rank+1], dtype=float) / comm.size

    if 0: #XXX which one is best?
        rel_j = np.cumsum(np.hstack((0,cost_j)))/np.sum(cost_j)
        ja,jb = np.argmin(np.abs(rel_j-a)), np.argmin(np.abs(rel_j-b))
    else:
        rel_j = np.cumsum(cost_j)/np.sum(cost_j)
        ja,jb = np.take(np.argwhere(np.bitwise_and(a<rel_j,rel_j<=b)),[0,-1])

    mypart = np.array([cost_j[ja:jb+1].sum()/np.sum(cost_j)])
    parts_r = np.empty(comm.size, dtype=float)
    comm.all_gather(mypart, parts_r)
    assert np.abs(np.sum(parts_r)-1) < 1e-9, (parts_r,np.sum(parts_r))
    comm.barrier()
    return iter(lm_j[ja:jb+1])

