from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from gpaw.mpi import serial_comm
from gpaw.test import equal
from gpaw.xc.rpa import RPACorrelation
from gpaw.xc.fxc import FXCCorrelation

a0 = 5.43
Ni = bulk('Ni', 'fcc')
Ni.set_initial_magnetic_moments([0.7])

calc = GPAW(mode='pw',
            kpts=(3, 3, 3),
            occupations=FermiDirac(0.001),
            setups={'Ni': '10'},
            communicator=serial_comm)
Ni.set_calculator(calc)
E = Ni.get_potential_energy()
calc.diagonalize_full_hamiltonian(nbands=50)

rpa = RPACorrelation(calc, nfrequencies=8, skip_gamma=True)
E_rpa = rpa.calculate(ecut=[50])

fxc = FXCCorrelation(calc, nlambda=16, nfrequencies=8, skip_gamma=True)
E_fxc = fxc.calculate(ecut=[50])

equal(E_rpa, -7.826, 0.01)
equal(E_fxc, -7.826, 0.01)
