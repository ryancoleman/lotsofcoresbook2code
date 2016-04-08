import numpy as np
from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from gpaw.wavefunctions.pw import PW
from gpaw.response.gw import GW

a = 5.431
atoms = bulk('Si', 'diamond', a=a)

calc = GPAW(
            mode=PW(200),
            kpts=(3,3,3),
            xc='LDA',
            eigensolver='cg',
            occupations=FermiDirac(0.001),
            parallel={'band': 1},
            txt='Si_groundstates.txt'
           )

atoms.set_calculator(calc)
atoms.get_potential_energy()

calc.diagonalize_full_hamiltonian()
calc.write('Si_groundstate.gpw','all')

gw = GW(
        file='Si_groundstate.gpw',
        nbands=None,
        bands=np.array([2,3,4,5]),
        kpoints=None,
        ecut=100.,
        txt='Si_EXX.out'
       )

gw.get_exact_exchange()

for wlin in [25.,50.,75.,100.]:

    for dw in [0.02,0.05,0.1,0.2,0.5]:

        gw = GW(
                file='Si_groundstate.gpw',
                nbands=None,
                bands=np.array([2,3,4,5]),
                kpoints=None,
                ecut=100.,
                w=np.array([wlin, 150., dw]),
                wpar=4,
                txt='Si_GW_wlin%s_dw%s.out' % (wlin, dw)
               )

        gw.get_QP_spectrum(file='Si_GW_wlin%s_dw%s.pckl' % (wlin, dw))
