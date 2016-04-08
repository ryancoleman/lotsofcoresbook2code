import numpy as np
from gpaw.mpi import size, rank, world
from gpaw.response.parallel import par_write, par_read

Nw = 400
npw = 10
assert Nw % size == 0
Nw_local = Nw // size

chi0_wGG = np.ones((Nw_local, npw, npw),dtype=complex) * rank

filename = 'chi0'
name = 'chi0_wGG'
par_write(filename, name, world, chi0_wGG)

chi0_wGG_new = par_read(filename, name)

assert (np.abs(chi0_wGG - chi0_wGG_new) < 1e-10).all()
