from os import system
from ase import *
from ase.lattice import bulk
from ase.dft.kpoints import monkhorst_pack
from gpaw import *
from gpaw.test import equal
from gpaw.xc.fxc import FXCCorrelation
from gpaw.mpi import world, serial_comm

if world.rank == 0:
    a0  = 5.43
    Ni = bulk('Ni', 'fcc')
    Ni.set_initial_magnetic_moments([0.7])

    kpts = monkhorst_pack((3,3,3))

    calc = GPAW(mode='pw',
                kpts=kpts,
                occupations=FermiDirac(0.001),
                setups={'Ni': '10'},
                communicator=serial_comm)
    
    Ni.set_calculator(calc)
    Ni.get_potential_energy()
    calc.diagonalize_full_hamiltonian()
    calc.write('Ni.gpw', mode='all')

world.barrier()

rpa = FXCCorrelation('Ni.gpw', xc='RPA',
                     nfrequencies=8, skip_gamma=True)
E_rpa = rpa.calculate(ecut=[50])

ralda = FXCCorrelation('Ni.gpw', xc='rALDA', unit_cells=[2,1,1],
                       nfrequencies=8, skip_gamma=True)
E_ralda = ralda.calculate(ecut=[50])

rapbe = FXCCorrelation('Ni.gpw', xc='rAPBE', unit_cells=[2,1,1],
                       nfrequencies=8, skip_gamma=True)
E_rapbe = rapbe.calculate(ecut=[50])

equal(E_rpa, -7.827, 0.01)
equal(E_ralda, -7.501, 0.01)
equal(E_rapbe, -7.444, 0.01)
