from __future__ import print_function
from ase.optimize import BFGS
from ase.structure import molecule
from ase.parallel import paropen
from gpaw import GPAW, PW
from gpaw.xc.exx import EXX

# N -------------------------------------------

N = molecule('N')
N.cell = (6, 6, 7)
calc = GPAW(mode=PW(600),
            dtype=complex,
            nbands=8,
            maxiter=300,
            xc='PBE',
            hund=True,
            txt='N_pbe.txt',
            convergence={'density': 1.e-6})

N.calc = calc
E1_pbe = N.get_potential_energy()

exx = EXX(calc, txt='N_exx.txt')
exx.calculate()
E1_hf = exx.get_total_energy()

calc.diagonalize_full_hamiltonian(nbands=4800)
calc.write('N.gpw', mode='all')

# N2 ------------------------------------------

N2 = molecule('N2')
N2.cell = (6, 6, 7)
calc = GPAW(mode=PW(600),
            dtype=complex,
            maxiter=300,
            xc='PBE',
            txt='N2_pbe.txt',
            convergence={'density': 1.e-6})

N2.calc = calc
dyn = BFGS(N2)
dyn.run(fmax=0.05)
E2_pbe = N2.get_potential_energy()

exx = EXX(calc, txt='N2_exx.txt')
exx.calculate()
E2_hf = exx.get_total_energy()

f = paropen('PBE_HF.dat', 'w')
print('PBE: ', E2_pbe - 2 * E1_pbe, file=f)
print('HF: ', E2_hf - 2 * E1_hf, file=f)
f.close()

calc.diagonalize_full_hamiltonian(nbands=4800)
calc.write('N2.gpw', mode='all')
