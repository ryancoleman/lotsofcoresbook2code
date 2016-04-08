from __future__ import print_function
# Refer to J. Harl, et.al, PRB 81, 115126 (2010)
#   Table III for parameter used in this calculation
#   Fig. 2 to compare the correlation energy.

import numpy as np
from gpaw import GPAW
from gpaw.xc.rpa_correlation_energy import RPACorrelation
from gpaw.mpi import serial_comm, size
from ase.parallel import paropen

tag = 'Nabulk'
f = paropen('%s_RPA.dat' %(tag), 'a')

for i in range(9):
    ecut = 120 + 10 * i

    calc = GPAW('gs_%s.gpw' %(tag), communicator=serial_comm, txt=None)
    rpa = RPACorrelation(calc, txt='rpa_%s_%s.txt' %(tag,ecut))
    E = rpa.get_rpa_correlation_energy(ecut=ecut,
                                       skip_gamma=False,
                                       directions=[[0,1.0]],
                                       kcommsize=size,
                                       dfcommsize=size)

    print('%s  '%(ecut), E, file=f)
    
f.close()
