from ase import *
from ase.structure import molecule
from gpaw import *
from gpaw.mpi import serial_comm
from gpaw.test import equal
from gpaw.xc.rpa_correlation_energy import RPACorrelation

ecut = 50

He = Atoms('He')
He.center(vacuum=1.0)

calc = GPAW(mode='pw',
            dtype=complex,
            xc='PBE',
            communicator=serial_comm)
He.set_calculator(calc)
He.get_potential_energy()
calc.diagonalize_full_hamiltonian()

rpa = RPACorrelation(calc)
C6_rpa, C6_0 = rpa.get_C6_coefficient(ecut=ecut,
                                      direction=2)

equal(C6_0, 1.772, 0.01)
equal(C6_rpa, 1.387, 0.01)
