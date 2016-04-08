from __future__ import print_function
import numpy as np
from gpaw import GPAW
from gpaw.dipole_correction import DipoleCorrection
from gpaw.poisson import PoissonSolver
from ase.structure import molecule
from gpaw.mpi import rank
from gpaw.utilities import h2gpts


"""
Test the dipole correction code by comparing this system:

    H
z1   O  z2
    H


(where z1 and z2 denote points where the potential is probed)

Expected potential:

      -----
     /
    /
----

to this system:

   H           H
z1  O    z2   O
   H           H


Expected potential:

       -------
      /       \
     /         \
-----           ------

The height of the two potentials are tested to be the same.

Enable if-statement in the bottom for nice plots
"""


system1 = molecule('H2O')
system1.set_pbc((True, True, False))
system1.cell = 4.0 * np.array([[1.0, -1.5, 0.0], [1.0, 1.0, 0.0],
                               [0., 0., 1.]])
system1.center(vacuum=10.0, axis=2)

system2 = system1.copy()
system2.positions *= [1.0, 1.0, -1.0]
system2 += system1
system2.center(vacuum=6.0, axis=2)

convergence = dict(density=1e-5)

calc1 = GPAW(mode='lcao',
             convergence=convergence,
             gpts=h2gpts(0.25, system1.cell, idiv=8),
             poissonsolver=DipoleCorrection(PoissonSolver(relax='GS',
                                                          eps=1e-11), 2))

system1.set_calculator(calc1)
system1.get_potential_energy()
v1 = calc1.get_effective_potential(pad=False)



calc2 = GPAW(mode='lcao',
             convergence=convergence,
             gpts=h2gpts(0.25, system2.cell, idiv=8),
             poissonsolver=PoissonSolver(relax='GS', eps=1e-11))

system2.set_calculator(calc2)
system2.get_potential_energy()
v2 = calc2.get_effective_potential(pad=False)


def get_avg(v):
    nx, ny, nz = v.shape
    vyz = v.sum(axis=0) / nx
    vz = vyz.sum(axis=0) / ny
    return vz, vyz

if rank == 0:
    vz1, vyz1 = get_avg(v1)
    vz2, vyz2 = get_avg(v2)

    # Compare values that are not right at the end of the array
    # (at the end of the array things can "oscillate" a bit)
    dvz1 = vz1[-5] - vz1[4]
    dvz2 = vz2[4] - vz2[len(vz2) // 2]
    print(dvz1, dvz2)
    
    err1 = abs(dvz1 - dvz2)
    
    # Comparison to what the values were when this test was last modified:
    ref_value = 2.07342988218
    err2 = abs(dvz1 - ref_value)

    if 0:
        import pylab as pl
        pl.imshow(vyz1)
        pl.figure()
        pl.imshow(vyz2)
        pl.figure()
        pl.plot(vz1)
        pl.plot(vz2)
        pl.show()

    print('Ref value of previous calculation', ref_value)
    print('Value in this calculation', dvz1)

    # fine grid needed to achieve convergence!
    print('Error', err1, err2)
    assert err1 < 5e-3, err1
    assert err2 < 2e-4, err2
