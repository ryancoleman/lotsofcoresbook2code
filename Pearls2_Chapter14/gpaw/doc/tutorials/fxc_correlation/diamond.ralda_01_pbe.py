from __future__ import print_function
from ase import Atoms
from ase.lattice import bulk
from ase.dft import monkhorst_pack
from ase.parallel import paropen
from gpaw import GPAW, FermiDirac
from gpaw.wavefunctions.pw import PW
from gpaw.xc.exx import EXX
import numpy as np

# Monkhorst-Pack grid shifted to be gamma centered
k = 8
kpts = monkhorst_pack([k, k, k])
kpts += [1. / (2 * k), 1. / (2 * k), 1. / (2 * k)]

cell = bulk('C', 'fcc', a=3.553).get_cell()
a = Atoms('C2', cell=cell, pbc=True,
          scaled_positions=((0, 0, 0), (0.25, 0.25, 0.25)))

calc = GPAW(mode=PW(600),
            xc='PBE',
            occupations=FermiDirac(width=0.01),
            convergence={'density': 1.e-6},
            kpts=kpts,
            txt='diamond.ralda_01_pbe.txt',
            )

a.set_calculator(calc)
E_pbe = a.get_potential_energy()

exx = EXX(calc, txt='diamond.ralda_01_exx.txt')
exx.calculate()
E_hf = exx.get_total_energy()

E_C = np.loadtxt('CO.ralda.PBE_HF_C.dat')

f = paropen('diamond.ralda.PBE_HF_diamond.dat', 'w')
print('PBE: ', E_pbe / 2 - E_C[0], file=f)
print('HF: ', E_hf / 2 - E_C[1], file=f)
f.close()

calc.diagonalize_full_hamiltonian()
calc.write('diamond.ralda.pbe_wfcs.gpw', mode='all')
