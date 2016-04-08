from __future__ import print_function
import numpy as np
from time import time, ctime
from datetime import timedelta
from ase.lattice import bulk
from ase.units import Hartree
from gpaw import GPAW, FermiDirac
from gpaw.response.gw import GW

starttime = time()

a = 3.567
atoms = bulk('C', 'diamond', a=a)

kpts = (2,2,2)

calc = GPAW(
            h=0.20,
            kpts=kpts,
            xc='LDA',
            txt='C_gs.txt',
            nbands=30,
            symmetry='off',
            convergence={'bands':25},
            eigensolver='cg',
            occupations=FermiDirac(0.001)
           )

atoms.set_calculator(calc)
atoms.get_potential_energy()
calc.write('C_gs.gpw','all')

file='C_gs.gpw'

gw = GW(
        file=file,
        nbands=25,
        bands=np.array([3,4]),
        w=np.array([75., 100., 0.05]),
        ecut=25.,
        eta=0.2,
        ppa=False,
        hilbert_trans=False,
       )

gw.get_exact_exchange()
gw.get_QP_spectrum()

QP_False = gw.QP_skn * Hartree

gw = GW(
        file=file,
        nbands=25,
        bands=np.array([3,4]),
        w=np.array([75., 100., 0.05]),
        ecut=25.,
        eta=0.2,
        ppa=True,
        hilbert_trans=False,
       )

gw.get_QP_spectrum()

QP_True = gw.QP_skn * Hartree

if not (np.abs(QP_False - QP_True) < 0.05).all():
    raise AssertionError("dynamic GW not equal to PPA")

totaltime = round(time() - starttime)
print("GW test finished in %s " %(timedelta(seconds=totaltime)))
