from ase import Atoms
from gpaw import GPAW, FermiDirac
from gpaw.eigensolvers.rmm_diis_old import RMM_DIIS
from gpaw.mixer import MixerSum
from gpaw.test import equal

a = 2.87
bulk = Atoms('Fe2',
             scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)],
             magmoms=[2.20, 2.20],
             cell=(a, a, a),
             pbc=True)
mom0 = sum(bulk.get_initial_magnetic_moments())
h = 0.2
conv = {'eigenstates': 0.1, 'density':0.1, 'energy':0.01}
calc = GPAW(h=h,
            eigensolver=RMM_DIIS(),
            mixer=MixerSum(0.1,3),
            nbands=11,
            kpts=(3, 3, 3),
            convergence=conv,
            occupations=FermiDirac(0.1, fixmagmom=True))
bulk.set_calculator(calc)
e = bulk.get_potential_energy()
niter = calc.get_number_of_iterations()
mom = calc.get_magnetic_moment()
equal(mom, mom0, 1e-5)

energy_tolerance = 0.0002
niter_tolerance = 0
equal(e, -20.3251, energy_tolerance)
