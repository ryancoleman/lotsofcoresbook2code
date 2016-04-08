from __future__ import print_function
from ase.structure import molecule
from gpaw import GPAW
from gpaw import dscf
from gpaw.test import equal

# Ground state calculation
#------------------------------------------------------------------

calc_mol = GPAW(nbands=8, h=0.2, xc='PBE', spinpol=True,
                convergence={'energy': 100,
                             'density': 100,
                             #'eigenstates': 1.0e-9,
                             'bands': -1})

CO = molecule('CO')
CO.center(vacuum=3)
CO.set_calculator(calc_mol)
E_gs = CO.get_potential_energy()
niter_gs = calc_mol.get_number_of_iterations()

'''Get the pseudowavefunctions and projector overlaps of the
   state which is to be occupied. n=5,6 is the 2pix and 2piy orbitals'''
n = 5
molecule = [0,1]
wf_u = [kpt.psit_nG[n] for kpt in calc_mol.wfs.kpt_u]
p_uai = [dict([(molecule[a], P_ni[n]) for a, P_ni in kpt.P_ani.items()])
         for kpt in calc_mol.wfs.kpt_u]

# Excited state calculations
#--------------------------------------------
calc_1 = GPAW(nbands=8, h=0.2, xc='PBE', spinpol=True,
              convergence={'energy': 100,
                           'density': 100,
                           #'eigenstates': 1.0e-9,
                           'bands': -1})
CO.set_calculator(calc_1)
weights = {0: [0.,0.,0.,1.], 1: [0.,0.,0.,-1.]}
lumo = dscf.MolecularOrbital(calc_1, weights=weights)
dscf.dscf_calculation(calc_1, [[1.0, lumo, 1]], CO)
E_es1 = CO.get_potential_energy()
niter_es1 = calc_1.get_number_of_iterations()
calc_1.write('dscf_CO_es1.gpw', mode='all')

calc_2 = GPAW(nbands=8, h=0.2, xc='PBE', spinpol=True,
              convergence={'energy': 100,
                          'density': 100,
                           #'eigenstates': 1.0e-9,
                           'bands': -1})
CO.set_calculator(calc_2)
lumo = dscf.AEOrbital(calc_2, wf_u, p_uai)
dscf.dscf_calculation(calc_2, [[1.0, lumo, 1]], CO)
E_es2 = CO.get_potential_energy()
niter_es2 = calc_2.get_number_of_iterations()
calc_2.write('dscf_CO_es2.gpw', mode='all')

equal(E_es1, E_gs+5.8, 0.1)
equal(E_es1, E_es2, 0.001)

energy_tolerance = 0.001
niter_tolerance = 2
print(E_gs, niter_gs)
print(E_es1, niter_es1)
print(E_es2, niter_es2)
equal(E_gs, -15.1924620949, energy_tolerance)
equal(E_es1, -9.36671359062, energy_tolerance)
equal(E_es2, -9.3667912622, energy_tolerance)
