import os
from gpaw import GPAW, restart
from ase import Atoms
from gpaw.test import equal
from gpaw.mpi import world, rank
from math import sqrt
import numpy as np

# Test the reading of wave functions as file references

modes = ['gpw']
try:
    import _gpaw_hdf5
    modes.append('hdf5')
except ImportError:
    pass

d = 3.0
atoms = Atoms('Na3', positions=[( 0, 0, 0),
                              ( 0, 0, d),
                              ( 0, d*sqrt(3./4.), d/2.)],
                   magmoms=[1.0, 1.0, 1.0],
                   cell=(3.5, 3.5, 4.+2/3.),
                   pbc=True)

# Only a short, non-converged calcuation
conv = {'eigenstates': 1.24, 'energy':2e-1, 'density':1e-1}
calc = GPAW(h=0.30, kpts=(1,1,3), 
            setups={'Na': '1'},
            nbands=3, convergence=conv)
atoms.set_calculator(calc)
e0 = atoms.get_potential_energy()
wf0 = calc.get_pseudo_wave_function(2, 1, 1, broadcast=True)
# Write the restart file(s)
for mode in modes:
    calc.write('tmp.%s' % mode, 'all')

del calc
# Now read with single process
comm = world.new_communicator(np.array((rank,)))
for mode in modes:
    calc = GPAW('tmp.%s' % mode, communicator=comm)
    wf1 = calc.get_pseudo_wave_function(2, 1, 1)
    diff = np.abs(wf0 - wf1)
    assert(np.all(diff < 1e-12))

world.barrier()
if rank == 0:
    for mode in modes:
        os.remove('tmp.%s' % mode)
