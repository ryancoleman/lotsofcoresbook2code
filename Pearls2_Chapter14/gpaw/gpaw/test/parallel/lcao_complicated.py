from __future__ import print_function
import numpy as np
from ase.lattice.surface import fcc111

from gpaw.mpi import world
from gpaw import GPAW, PoissonSolver
from gpaw.utilities import compiled_with_sl

# This test verifies that energy and forces are (approximately)
# parallelization independent
#
# Tests the LCAO energy and forces in non-orthogonal cell with
# simultaneous parallelization over bands, domains and k-points (if
# enough CPUs are available), both with and without scalapack
# (if scalapack is available).
#
# Run with 1, 2, 4 or 8 (best) CPUs.
#
# This test covers many cases not caught by lcao_parallel or
# lcao_parallel_kpt
#
# Written November 24, 2011, r8567

system = fcc111('Au', size=(1, 3, 1))
system.numbers[0] = 8
# It is important that the number of atoms is uneven; this
# tests the case where the band parallelization does not match
# the partitioning of orbitals between atoms (the middle atom has orbitals
# on distinct band descriptor ranks)

system.center(vacuum=3.5, axis=2)
system.rattle(stdev=0.2, seed=17)
#from ase.visualize import view
#view(system)

#system.set_pbc(0)
#system.center(vacuum=3.5)
from gpaw import FermiDirac

def calculate(parallel, comm=world, Eref=None, Fref=None):
    calc = GPAW(mode='lcao',
                basis=dict(O='dzp', Au='sz(dzp)'),
                occupations=FermiDirac(0.1),
                kpts=(4, 1, 1),
                #txt=None,
                communicator=comm,
                poissonsolver=PoissonSolver(relax='J', eps=1e-8),
                nbands=16,
                parallel=parallel,
                h=0.35)
    system.set_calculator(calc)
    E = system.get_potential_energy()
    F = system.get_forces()

    if world.rank == 0:
        print('Results')
        print('-----------')
        print(E)
        print(F)
        print('-----------')

    if Eref is not None:
        Eerr = abs(E - Eref)
        assert Eerr < 1e-8, 'Bad E: err=%f; parallel=%s' % (Eerr, parallel)
    if Fref is not None:
        Ferr = np.abs(F - Fref).max()
        assert Ferr < 1e-6, 'Bad F: err=%f; parallel=%s' % (Ferr, parallel)
    return E, F

# First calculate reference energy and forces E and F
#
# If we want to really dumb things down, enable this to force an
# entirely serial calculation:
if 0: 
    serial = world.new_communicator([0])
    E = 0.0
    F = np.zeros((len(system), 3))
    if world.rank == 0:
        E, F = calculate({}, serial)
    E = world.sum(E)
    world.sum(F)
else:
    # Normally we'll just do it in parallel;
    # that case is covered well by other tests, so we can probably trust it
    E, F = calculate({}, world)

def check(parallel):
    return calculate(parallel, comm=world, Eref=E, Fref=F)

assert world.size in [1, 2, 4, 8], ('Number of CPUs %d not supported'
                                    % world.size)

parallel = dict(domain=1, band=1)
sl_cpus = world.size
if world.size % 2 == 0:
    parallel['band'] = 2
if world.size % 4 == 0:
    parallel['domain'] = 2

# If size is 8, this will also use kpt parallelization.  This test should
# run with 8 CPUs for best coverage of parallelizations
if world.size == 8:
    sl_cpus = 4

if world.size > 1:
    check(parallel)

if compiled_with_sl():
    parallel['sl_auto'] = True
    check(parallel)
