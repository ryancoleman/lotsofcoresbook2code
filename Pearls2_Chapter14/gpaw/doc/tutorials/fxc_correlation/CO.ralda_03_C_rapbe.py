from __future__ import print_function
from ase.parallel import paropen
from ase.units import Hartree
from gpaw.xc.rpa import RPACorrelation
from gpaw.xc.fxc import FXCCorrelation
from gpaw.mpi import world

fxc0 = FXCCorrelation('CO.ralda.pbe_wfcs_C.gpw',
                      xc='rAPBE',
                      txt='CO.ralda_03_C_rapbe.txt',
                      wcomm=world.size)

E0_i = fxc0.calculate(ecut=400)

f = paropen('CO.ralda_rapbe_C.dat', 'w')
for ecut, E0 in zip(fxc0.ecut_i, E0_i):
    print(ecut * Hartree, E0, file=f)
f.close()

rpa0 = RPACorrelation('CO.ralda.pbe_wfcs_C.gpw',
                      txt='CO.ralda_03_C_rpa.txt',
                      wcomm=world.size)

E0_i = rpa0.calculate(ecut=400)

f = paropen('CO.ralda_rpa_C.dat', 'w')
for ecut, E0 in zip(rpa0.ecut_i, E0_i):
    print(ecut * Hartree, E0, file=f)
f.close()
