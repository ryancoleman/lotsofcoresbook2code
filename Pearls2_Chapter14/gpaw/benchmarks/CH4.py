from ase.structure import molecule
from gpaw import GPAW
from gpaw.eigensolvers import RMM_DIIS
from gpaw.mpi import rank

atoms = molecule('CH4')
atoms.set_pbc(True)
atoms.center(3.5)
calc = GPAW(h=0.2, convergence={'eigenstates' : 1e-5},
            eigensolver=RMM_DIIS(keep_htpsit=True))#, txt='out.txt')
atoms.set_calculator(calc)
e = atoms.get_potential_energy()
