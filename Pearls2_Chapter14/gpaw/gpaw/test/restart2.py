"""Test automatically write out of restart files"""

import os
from glob import glob
from gpaw import GPAW
from ase import Atom, Atoms
from gpaw.test import equal

restart = 'gpaw-restart'
result  = 'gpaw-result'

# H atom:
H = Atoms([Atom('H')])
H.center(vacuum=3.0)

calc = GPAW(gpts=(32, 32, 32), nbands=1)
calc.attach(calc.write, 4, restart)
H.set_calculator(calc)
e = H.get_potential_energy()
niter = calc.get_number_of_iterations()
calc.write(result)



# the two files should be equal
from gpaw.mpi import rank
if rank == 0:
    for f in ['gpaw-restart', 'gpaw-result']:
        os.system('rm -rf %s; mkdir %s; cd %s; tar xf ../%s.gpw' %
                  (f, f, f, f))

    # make a list of the files to compare
    files_restart = glob(restart+'/*')
    files_result = glob(result+'/*')

    for f1, f2 in zip(files_restart, files_result):
        # cmp is a byte-by-byte comparison, more portable than diff
        assert os.system('cmp %s %s' % (f1, f2)) == 0
    os.system('rm -rf gpaw-restart gpaw-result')

energy_tolerance = 0.00006
niter_tolerance = 0
equal(e, 0.0451789, energy_tolerance)
