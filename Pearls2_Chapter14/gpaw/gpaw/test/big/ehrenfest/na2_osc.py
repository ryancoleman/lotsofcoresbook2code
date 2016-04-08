from __future__ import print_function

import os, time
import numpy as np

from ase import Atoms
from ase.parallel import paropen
from ase.units import Hartree, Bohr
from ase.io.trajectory import PickleTrajectory
from ase.calculators.singlepoint import SinglePointCalculator
from gpaw import GPAW
from gpaw.mpi import world
from gpaw.tddft import TDDFT
from gpaw.tddft.ehrenfest import EhrenfestVelocityVerlet

# -------------------------------------------------------------------

name = 'na2_osc'

# Equilibrium distance in Ang cf. setups page for Na dimer
d_bond = 3.29 # ~3.0 during oscillation
d_disp = 0.1

# Timestep and expected oscillatory period in attoseconds
timestep = 5.0
period = 2.1e5 # ~19.7 meV cf. CRC Handbook of Phys. & Chem. #09_08_91

ndiv = int(np.ceil(0.1e3 / timestep)) # update stats every 0.1 fs
niter = ndiv * int(np.ceil(period / (ndiv * timestep)))

if __name__ == '__main__':
    if not os.path.isfile(name + '_gs.gpw'):
        atoms = Atoms('Na2', positions=[(0, 0, 0), (0, 0, d_bond + d_disp)])
        atoms.set_pbc(False)
        atoms.center(vacuum=6.0)
        cell_c = np.sum(atoms.get_cell()**2, axis=1)**0.5
        N_c = 16 * np.round(cell_c / (0.25 * 16))
        calc = GPAW(gpts=N_c, nbands=1, basis='dzp', setups={'Na': '1'},
                    txt=name + '_gs.txt')
        atoms.set_calculator(calc)
        atoms.get_potential_energy()
        calc.write(name + '_gs.gpw', mode='all')
        del atoms, calc
        time.sleep(10)

    while not os.path.isfile(name + '_gs.gpw'):
        print('Node %d waiting for file...' % world.rank)
        time.sleep(10)
    world.barrier()

    tdcalc = TDDFT(name + '_gs.gpw', txt=name + '_td.txt', propagator='EFSICN')
    ehrenfest = EhrenfestVelocityVerlet(tdcalc)
    traj = PickleTrajectory(name + '_td.traj', 'w', tdcalc.get_atoms())

    t0 = time.time()
    f = paropen(name + '_td.log', 'w')
    for i in range(1, niter+1):
        ehrenfest.propagate(timestep)

        if i % ndiv == 0:
            rate = 60 * ndiv / (time.time()-t0)
            ekin = tdcalc.atoms.get_kinetic_energy()
            epot = tdcalc.get_td_energy() * Hartree
            F_av = ehrenfest.F * Hartree / Bohr
            print('i=%06d (%6.2f min^-1), ekin=%13.9f, epot=%13.9f, etot=%13.9f' % (i, rate, ekin, epot, ekin+epot), file=f)
            t0 = time.time()

            # Hack to prevent calls to GPAW::get_potential_energy when saving
            spa = tdcalc.get_atoms()
            spc = SinglePointCalculator(epot, F_av, None, None, spa)
            spa.set_calculator(spc)
            traj.write(spa)
    f.close()
    traj.close()
