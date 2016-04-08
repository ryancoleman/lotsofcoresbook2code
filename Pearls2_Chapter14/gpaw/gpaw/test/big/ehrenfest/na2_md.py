from __future__ import print_function

import os, time
import numpy as np

from ase import Atoms
from ase.parallel import paropen
from ase.units import fs, Hartree, Bohr
from ase.io.trajectory import PickleTrajectory
from ase.md.verlet import VelocityVerlet
from gpaw import GPAW
from gpaw.mpi import world

# -------------------------------------------------------------------

name = 'na2_md'

# Equilibrium distance in Ang cf. setups page for Na dimer
d_bond = 3.29 # ~3.0 during oscillation
d_disp = 0.1

# Timestep and expected oscillatory period in attoseconds
timestep = 500.0
period = 2.1e5 # ~19.7 meV cf. CRC Handbook of Phys. & Chem. #09_08_91

ndiv = int(np.ceil(10e3 / timestep)) # update stats every 10 fs
niter = ndiv * int(np.ceil(3 * period / (ndiv * timestep)))

class Timing:
    def __init__(self, f):
        self.f = f
        self.t0 = time.time()
        self.i = 0

    def __call__(self, atoms):
        rate = 60 * ndiv / (time.time()-self.t0)
        ekin = atoms.get_kinetic_energy()
        epot = atoms.get_potential_energy()
        print('i=%06d (%6.2f min^-1), ekin=%13.9f, epot=%13.9f, etot=%13.9f' % (self.i, rate, ekin, epot, ekin+epot), file=self.f)
        self.t0 = time.time()
        self.i += ndiv

if __name__ == '__main__':
    if not os.path.isfile(name + '_gs.gpw'):
        atoms = Atoms('Na2', positions=[(0, 0, 0), (0, 0, d_bond + d_disp)])
        atoms.set_pbc(False)
        atoms.center(vacuum=6.0)
        atoms.set_velocities(np.zeros_like(atoms.get_positions()))
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

    tdcalc = GPAW(name + '_gs.gpw', txt=name + '_td.txt')
    tdcalc.forces.reset() #XXX debug
    tdcalc.initialize_positions()
    atoms = tdcalc.get_atoms()

    traj = PickleTrajectory(name + '_td.traj', 'w', atoms)
    verlet = VelocityVerlet(atoms, timestep * 1e-3 * fs,
                            logfile=paropen(name + '_td.verlet', 'w'),
                            trajectory=traj)
    verlet.attach(Timing(paropen(name + '_td.log', 'w')), ndiv, atoms)
    verlet.run(niter)
    traj.close()
