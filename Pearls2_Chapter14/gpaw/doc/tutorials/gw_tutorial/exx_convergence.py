import numpy as np
from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from gpaw.wavefunctions.pw import PW
from gpaw.response.gw import GW

a = 5.431
atoms = bulk('Si', 'diamond', a=a)

for k in [3,5,7,9]:

    kpts = (k,k,k)

    for ecut in [50,100,150,200]:

        calc = GPAW(
                    mode=PW(ecut),
                    kpts=kpts,
                    xc='LDA',
                    eigensolver='cg',
                    occupations=FermiDirac(0.001),
                    parallel={'band': 1},
                    txt='Si_groundstate_k%s_ecut%s.txt' % (k, ecut)
                   )

        atoms.set_calculator(calc)
        atoms.get_potential_energy()

        calc.diagonalize_full_hamiltonian()
        calc.write('Si_groundstate_k%s_ecut%s.gpw' % (k, ecut),'all')

        gw = GW(
                file='Si_groundstate_k%s_ecut%s.gpw' % (k, ecut),
                nbands=None,
                bands=np.array([2,3,4,5]),
                kpoints=None,
                ecut=ecut,
                ppa=True,
                txt='Si_EXX_k%s_ecut%s.out' % (k, ecut)
               )

        gw.get_exact_exchange(file='Si_EXX_k%s_ecut%s.pckl' % (k, ecut))
