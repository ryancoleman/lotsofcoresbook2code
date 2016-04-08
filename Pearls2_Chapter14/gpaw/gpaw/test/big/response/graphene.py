from __future__ import print_function
import numpy as np
import sys
import time

from math import sqrt
from ase import Atoms, Atom, view
from ase.units import Bohr
from ase.parallel import paropen
from gpaw import GPAW
from gpaw.response.grid_chi import CHI
from gpaw.mpi import serial_comm, rank, size
from gpaw.utilities import devnull


if rank != 0:
  sys.stdout = devnull 

GS1 = 1
GS2 = 1
EELS = 1
OpticalLimit = 1

nband = 40

if GS1:

    kpts = (30,30,1)
    a=1.42

    atoms = Atoms([Atom('C',(1/3.0,1/3.0,0)),
                   Atom('C',(2/3.0,2/3.0,0))],
                   pbc=(1,1,1))
    atoms.set_cell([(sqrt(3)*a/2.0,3/2.0*a,0),
                    (-sqrt(3)*a/2.0,3/2.0*a,0),
                    (0,0,30)],
                   scale_atoms=True)
    atoms.center(axis=2)
    calc = GPAW(h=0.2,
                xc='LDA',
                txt='out1.txt',
                kpts=kpts,
                basis='dzp',
                symmetry='off',
                nbands=nband+10,
                convergence={'bands':nband},
                eigensolver = 'cg',
                width=0.05)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    calc.write('graphene1.gpw','all')

if GS2:

    # Gs calculation shifted kpoint
    calc = GPAW('graphene1.gpw')
    q = np.array([0.00001, 0., 0.])
    kpts = q + calc.get_ibz_k_points()

    atoms = calc.atoms
    calc = GPAW(xc='LDA',
                kpts=kpts,
                h = 0.2,
                basis='dzp',
                symmetry='off',
                nbands=nband+10,
                convergence={'bands':nband},
                eigensolver = 'cg',
                width=0.05, txt = 'out2.txt')
    
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    calc.write('graphene2.gpw','all')   

if EELS:

    calc = GPAW('graphene1.gpw',communicator=serial_comm)
                
    dw = 0.1  
    wmax = 35.
    f = paropen('graphene_q_list', 'w')

    for i in range(1, 10):
        q = np.array([i/30., 0., 0.]) # Gamma - M
#        q = np.array([i/30., -i/30., 0.]) # Gamma - K 
        chi = CHI()
        chi.nband = nband
        chi.initialize((calc,), q, wmax, dw, eta=0.2, Ecut = 50 + (i-1)*10)
        chi.periodic()

        chi.get_EELS_spectrum('graphene_EELS_' + str(i))
        chi.check_sum_rule()

        print(sqrt(np.inner(chi.qq / Bohr, chi.qq / Bohr)), file=f)
    
if OpticalLimit:
    
    calc1 = GPAW('graphene1.gpw',communicator=serial_comm)
    calc2 = GPAW('graphene2.gpw',communicator=serial_comm)
                
    dw = 0.1  
    wmax = 35.

    q = np.array([0.00001, 0., 0.])
    chi = CHI()
    chi.nband = nband
    chi.initialize((calc1,calc2), q, wmax, dw, eta=0.2, Ecut = 40)
    chi.periodic()
    
    chi.get_absorption_spectrum('graphene_absorption')
    chi.check_sum_rule()
