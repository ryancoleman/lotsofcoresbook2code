from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from gpaw.mpi import serial_comm
from gpaw.test import equal
from gpaw.xc.rpa import RPACorrelation

a0 = 5.43
Si = bulk('Si', a=a0)

calc = GPAW(mode='pw',
            kpts={'size': (2, 2, 2), 'gamma': True},
            occupations=FermiDirac(0.001),
            communicator=serial_comm)
Si.set_calculator(calc)
E = Si.get_potential_energy()
calc.diagonalize_full_hamiltonian(nbands=50)

ecut = 50
rpa = RPACorrelation(calc, qsym=False, nfrequencies=8)
E_rpa_noqsym = rpa.calculate(ecut=[ecut])

rpa = RPACorrelation(calc, qsym=True, nfrequencies=8)
E_rpa_qsym = rpa.calculate(ecut=[ecut])

print(E_rpa_qsym, E_rpa_noqsym, 0.001)
equal(E_rpa_qsym, -12.61, 0.01)
