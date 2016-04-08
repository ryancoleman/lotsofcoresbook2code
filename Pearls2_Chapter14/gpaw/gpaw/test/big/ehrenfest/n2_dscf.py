from __future__ import print_function

import os, time
import numpy as np

from ase import Atoms
from ase.optimize.bfgs import BFGS
from ase.parallel import paropen
from ase.units import Hartree, Bohr
from ase.io.trajectory import PickleTrajectory
from ase.calculators.singlepoint import SinglePointCalculator
from gpaw import GPAW
from gpaw.mpi import world
from gpaw.occupations import FermiDirac
from gpaw.dscf import AEOrbital, dscf_calculation
from gpaw.utilities.dscftools import dscf_collapse_orbitals
from gpaw.tddft import TDDFT
from gpaw.tddft.ehrenfest import EhrenfestVelocityVerlet

name = 'n2_dscf'

# Equilibrium distance in Ang cf. setups page for N dimer
d_bond = 1.1
d_disp = 0.04 #XXX purely to get same cell as n2_osc

# Timestep and expected oscillatory period in attoseconds
timestep = 5.0
period = 1.924e4 # ~214.9 meV for $B^3 \Pi_g$ cf. Lofthus and Krupenie, 1977

ndiv = int(np.ceil(0.1e3 / timestep)) # update stats every 0.1 fs
niter = ndiv * int(np.ceil(2 * period / (ndiv * timestep)))

if __name__ == '__main__':
    if not os.path.isfile(name + '_gs.gpw'):
        # Calculator for quick ground state relaxation
        atoms = Atoms('N2', positions=[(0, 0, 0), (0, 0, d_bond + d_disp)])
        atoms.set_pbc(False)
        atoms.center(vacuum=6.0)
        cell_c = np.sum(atoms.get_cell()**2, axis=1)**0.5
        N_c = 8 * np.round(cell_c / (0.2 * 8))
        calc = GPAW(gpts=N_c, nbands=5, basis='dzp', #TODO xc='PBE'
                    txt=name + '_gs.txt', parallel={'band': 1})
        atoms.set_calculator(calc)

        # QuasiNewton relaxation before we increase the number of bands
        qn = BFGS(atoms, logfile=name + '_gs.log', trajectory=name + '_gs.traj')
        qn.run(0.01, steps=100)

        # Converge enough unoccupied bands for dSCF expansion
        calc.set(nbands=10, spinpol=True, occupations=FermiDirac(0.1),
                 convergence={'bands': -2, 'eigenstates': 1.0e-9})
        atoms.get_potential_energy()
        calc.write(name + '_gs.gpw', mode='all')
        del qn, calc, atoms
        time.sleep(10)

    while not os.path.isfile(name + '_gs.gpw'):
        print('Node %d waiting for %s...' % (world.rank, name + '_gs.gpw'))
        time.sleep(10)
    world.barrier()

    if not os.path.isfile(name + '_es.gpw'):
        # Resume ground state calculator and use as basis for excited state
        calc = GPAW(name + '_gs.gpw', txt=name + '_es.txt', parallel={'band': 1})
        calc.set_positions()
        atoms = calc.get_atoms()
        e1 = atoms.get_potential_energy()

        # Obtain the pseudowavefunctions and projector overlaps of the
        # state which is to be occupied. n=5,6 is the 2pix and 2piy orbitals
        n = 5
        molecule = [0, 1]
        wf_u = [kpt.psit_nG[n].copy() for kpt in calc.wfs.kpt_u]
        p_uai = [dict([(molecule[a], P_ni[n].copy()) for a, P_ni in 
                 kpt.P_ani.items()]) for kpt in calc.wfs.kpt_u]

        #calc = GPAW(h=0.2, nbands=10, xc='PBE', spinpol=True,
        #            convergence={'bands': -2, 'eigenstates': 1.0e-9},
        #            occupations=FermiDirac(0.1), txt=name + '_es.txt')
        #atoms.set_calculator(calc)

        lumo = AEOrbital(calc, wf_u, p_uai)
        dscf_calculation(calc, [[1.0, lumo, 1]], atoms)
        e2 = atoms.get_potential_energy()
        if world.rank == 0:
            print('e1:', e1, 'e2:', e2, 'de:', e2-e1)
        calc.write(name + '_es.gpw', mode='all')
        del wf_u, p_uai, lumo, calc, atoms
        time.sleep(10)

    while not os.path.isfile(name + '_es.gpw'):
        print('Node %d waiting for %s...' % (world.rank, name + '_es.gpw'))
        time.sleep(10)
    world.barrier()

    if not os.path.isfile(name + '_esx.gpw'):
        calc = GPAW(name + '_es.gpw', txt=name + '_esx.txt', parallel={'band': 1})
        calc.set_positions()
        dscf_collapse_orbitals(calc)
        calc.write(name + '_esx.gpw', mode='all')
        del calc
        time.sleep(10)

    while not os.path.isfile(name + '_esx.gpw'):
        print('Node %d waiting for %s...' % (world.rank, name + '_esx.gpw'))
        time.sleep(10)
    world.barrier()

    tdcalc = TDDFT(name + '_esx.gpw', txt=name + '_td.txt', propagator='EFSICN')
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
