from ase import Atoms
from gpaw import GPAW
from gpaw.wavefunctions.pw import PW
from gpaw.xc.fxc import FXCCorrelation
from gpaw.test import equal
from gpaw.mpi import world

if world.size == 1:
    scalapack1 = None
    scalapack2 = None
elif world.size == 2:
    scalapack1 = (2, world.size // 2, 32)    
    scalapack2 = None
else:
    scalapack1 = (2, world.size // 2, 32)
    scalapack2 = (2, world.size // 4, 32)

# H2 --------------------------------------
H2 = Atoms('H2', [(0,0,0),(0,0,0.7413)])
H2.set_pbc(True)
H2.set_cell((2., 2., 3.))
H2.center()
calc = GPAW(mode=PW(210),
            eigensolver='rmm-diis',
            dtype=complex,
            xc='LDA',
            basis='dzp',
            nbands=8,
            convergence={'density': 1.e-6})
H2.set_calculator(calc)
H2.get_potential_energy()
calc.diagonalize_full_hamiltonian(nbands=80, scalapack=scalapack1)
calc.write('H2.gpw', mode='all')

ralda = FXCCorrelation('H2.gpw', xc='rALDA')
E_ralda_H2 = ralda.calculate(ecut=[200])

rapbe = FXCCorrelation('H2.gpw', xc='rAPBE')
E_rapbe_H2 = rapbe.calculate(ecut=[200])

# H ---------------------------------------
H = Atoms('H', [(0,0,0)])
H.set_pbc(True)
H.set_cell((2., 2., 3.))
H.center()
calc = GPAW(mode=PW(210),
            eigensolver='rmm-diis',
            dtype=complex,
            xc='LDA',
            basis='dzp',
            nbands=4,
            hund=True,
            convergence={'density': 1.e-6})
H.set_calculator(calc)
H.get_potential_energy()
calc.diagonalize_full_hamiltonian(nbands=80, scalapack=scalapack2)
calc.write('H.gpw', mode='all')

ralda = FXCCorrelation('H.gpw', xc='rALDA')
E_ralda_H = ralda.calculate(ecut=[200])

rapbe = FXCCorrelation('H.gpw', xc='rAPBE')
E_rapbe_H = rapbe.calculate(ecut=[200])
                                      
# if rank == 0:
#    system('rm H2.gpw')
#    system('rm H.gpw')
#    system('rm fhxc_H2_rALDA_200_0.gpw')
#    system('rm fhxc_H_rALDA_200_0.gpw')
#    system('rm fhxc_H2_rAPBE_200_0.gpw')
#    system('rm fhxc_H_rAPBE_200_0.gpw')

equal(E_ralda_H2, -0.8411, 0.001)
equal(E_ralda_H, 0.002860, 0.00001)
equal(E_rapbe_H2, -0.7233, 0.001)
equal(E_rapbe_H, 0.016022, 0.00001)
