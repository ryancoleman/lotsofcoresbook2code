from ase.lattice import bulk
from gpaw import GPAW, FermiDirac, PW
from gpaw.mpi import serial_comm
from gpaw.xc.rpa import RPACorrelation
from gpaw.test import equal

bulk = bulk('Na', 'bcc', a=4.23)

ecut = 350
calc = GPAW(mode=PW(ecut), dtype=complex, basis='dzp',
            kpts={'size': (4, 4, 4), 'gamma': True},
            parallel={'domain': 1}, txt='gs_occ_pw.txt', nbands=4,
            occupations=FermiDirac(0.01),
            setups={'Na': '1'},
            )
bulk.set_calculator(calc)
bulk.get_potential_energy()
calc.write('gs_occ_pw.gpw')

calc = GPAW('gs_occ_pw.gpw', txt='gs_pw.txt', parallel={'band': 1})
calc.diagonalize_full_hamiltonian(nbands=520)
calc.write('gs_pw.gpw', 'all')

ecut = 120
calc = GPAW('gs_pw.gpw', communicator=serial_comm, txt=None)
rpa = RPACorrelation(calc, txt='rpa_%s.txt' % ecut)
E = rpa.calculate(ecut=[ecut])
equal(E, -1.106, 0.005)
