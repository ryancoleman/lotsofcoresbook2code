from __future__ import print_function
from ase import Atoms
from ase.parallel import paropen
from gpaw import GPAW, FermiDirac
from gpaw.mixer import MixerSum
from gpaw.xc.exx import EXX
from gpaw.wavefunctions.pw import PW

# CO ------------------------------------------

CO = Atoms('CO', [(0, 0, 0), (0, 0, 1.1283)])
CO.set_pbc(True)
CO.center(vacuum=3.0)
calc = GPAW(mode=PW(600),
            dtype=complex,
            xc='PBE',
            txt='CO.ralda_01_CO_pbe.txt',
            convergence={'density': 1.e-6})

CO.set_calculator(calc)
E0_pbe = CO.get_potential_energy()

exx = EXX(calc, txt='CO.ralda_01_CO_exx.txt')
exx.calculate()
E0_hf = exx.get_total_energy()

calc.diagonalize_full_hamiltonian()
calc.write('CO.ralda.pbe_wfcs_CO.gpw', mode='all')

# C -------------------------------------------

C = Atoms('C')
C.set_pbc(True)
C.set_cell(CO.cell)
C.center()
calc = GPAW(mode=PW(600),
            dtype=complex,
            xc='PBE',
            mixer=MixerSum(beta=0.1, nmaxold=5, weight=50.0),
            hund=True,
            occupations=FermiDirac(0.01, fixmagmom=True),
            txt='CO.ralda_01_C_pbe.txt',
            convergence={'density': 1.e-6})

C.set_calculator(calc)
E1_pbe = C.get_potential_energy()

exx = EXX(calc, txt='CO.ralda_01_C_exx.txt')
exx.calculate()
E1_hf = exx.get_total_energy()

f = paropen('CO.ralda.PBE_HF_C.dat', 'w')
print(E1_pbe, E1_hf, file=f)
f.close()

calc.diagonalize_full_hamiltonian()
calc.write('CO.ralda.pbe_wfcs_C.gpw', mode='all')

# O -------------------------------------------

O = Atoms('O')
O.set_pbc(True)
O.set_cell(CO.cell)
O.center()
calc = GPAW(mode=PW(600),
            dtype=complex,
            xc='PBE',
            mixer=MixerSum(beta=0.1, nmaxold=5, weight=50.0),
            hund=True,
            txt='CO.ralda_01_O_pbe.txt',
            convergence={'density': 1.e-6})

O.set_calculator(calc)
E2_pbe = O.get_potential_energy()

exx = EXX(calc, txt='CO.ralda_01_O_exx.txt')
exx.calculate()
E2_hf = exx.get_total_energy()

calc.diagonalize_full_hamiltonian()
calc.write('CO.ralda.pbe_wfcs_O.gpw', mode='all')

f = paropen('CO.ralda.PBE_HF_CO.dat', 'w')
print('PBE: ', E0_pbe - E1_pbe - E2_pbe, file=f)
print('HF: ', E0_hf - E1_hf - E2_hf, file=f)
f.close()
