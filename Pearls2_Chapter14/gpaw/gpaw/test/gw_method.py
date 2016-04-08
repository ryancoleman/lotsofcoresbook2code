from __future__ import print_function
import numpy as np
from time import time, ctime
from datetime import timedelta
from ase.lattice import bulk
from ase.units import Hartree
from gpaw import GPAW, FermiDirac
from gpaw.response.gw import GW
from gpaw.mpi import serial_comm

starttime = time()

a = 5.431
atoms = bulk('Si', 'diamond', a=a)

kpts = (2,2,2)

calc = GPAW(
            h=0.24,
            kpts=kpts,
            xc='LDA',
            txt='Si_gs.txt',
            nbands=10,
            convergence={'bands':8},
            occupations=FermiDirac(0.001)
           )

atoms.set_calculator(calc)
atoms.get_potential_energy()
calc.write('Si_gs.gpw','all')

file='Si_gs.gpw'
nbands=8
bands=np.array([3,4])
ecut=25./Hartree

gw = GW(
        file=file,
        nbands=8,
        bands=np.array([3,4]),
        w=np.array([10., 30., 0.05]),
        ecut=25.,
        eta=0.1,
        hilbert_trans=False
       )

gw.get_exact_exchange(communicator=serial_comm)

gw.get_QP_spectrum()

QP_False = gw.QP_skn * Hartree

gw = GW(
        file=file,
        nbands=8,
        bands=np.array([3,4]),
        w=np.array([10., 30., 0.05]),
        ecut=25.,
        eta=0.1,
        hilbert_trans=True
       )

gw.get_QP_spectrum()

QP_True = gw.QP_skn * Hartree

if not (np.abs(QP_False - QP_True) < 0.01).all():
    raise AssertionError("method 1 not equal to method 2")

totaltime = round(time() - starttime)
print("GW test finished in %s " %(timedelta(seconds=totaltime)))
