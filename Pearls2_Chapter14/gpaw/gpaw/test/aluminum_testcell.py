from __future__ import print_function
import numpy as np
import sys
import os
import time

from ase import Atom, Atoms
from ase.visualize import view
from ase.units import Bohr
from ase.lattice import bulk
from ase.utils import devnull

from gpaw import GPAW
from gpaw.eigensolvers.rmm_diis_old import RMM_DIIS
from gpaw.atom.basis import BasisMaker
from gpaw.response.df0 import DF
from gpaw.mixer import Mixer
from gpaw.mpi import serial_comm, rank, size


# Ground state calculation
a = 4.043
atoms = bulk('Al', 'fcc', a=a)
atoms.center()
calc = GPAW(gpts=(12,12,12),
            eigensolver=RMM_DIIS(),
            mixer=Mixer(0.1,3),
            kpts=(4,4,4),
            xc='LDA')

atoms.set_calculator(calc)
atoms.get_potential_energy()
calc.write('Al1.gpw','all')

# Excited state calculation
q = np.array([1./4.,0.,0.])
w = np.linspace(0, 24, 241)

df = DF(calc='Al1.gpw', q=q, w=w, eta=0.2, ecut=50)
#df.write('Al.pckl')
df.get_EELS_spectrum(filename='EELS_Al_1')

atoms = Atoms('Al8',scaled_positions=[(0,0,0),
                               (0.5,0,0),
                               (0,0.5,0),
                               (0,0,0.5),
                               (0.5,0.5,0),
                               (0.5,0,0.5),
                               (0.,0.5,0.5),
                               (0.5,0.5,0.5)],
              cell=[(0,a,a),(a,0,a),(a,a,0)],
              pbc=True)

calc = GPAW(gpts=(24,24,24),
            eigensolver=RMM_DIIS(),
            mixer=Mixer(0.1,3),
            kpts=(2,2,2),
            xc='LDA')

atoms.set_calculator(calc)
atoms.get_potential_energy()
calc.write('Al2.gpw','all')

# Excited state calculation
q = np.array([1./2.,0.,0.])
w = np.linspace(0, 24, 241)

df = DF(calc='Al2.gpw', q=q, w=w, eta=0.2, ecut=50)
#df.write('Al.pckl')
df.get_EELS_spectrum(filename='EELS_Al_2')

d1 = np.loadtxt('EELS_Al_1')
d2 = np.loadtxt('EELS_Al_2')
error1 = (d1[1:,1] - d2[1:,1]) / d1[1:,1] * 100
error2 = (d1[1:,2] - d2[1:,2]) / d1[1:,2] * 100

if error1.max() > 0.2 or error2.max() > 0.2: # percent
    print(error1.max(), error2.max())
    raise ValueError('Pls check spectrum !')

#if rank == 0:
#    os.remove('Al1.gpw')
#    os.remove('Al2.gpw')

