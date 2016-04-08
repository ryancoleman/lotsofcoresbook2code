# Note: Usesymm=None should be used for optical absorption spectra

import numpy as np
import sys

from math import sqrt
from ase import Atoms, Atom
from ase.units import Bohr
from ase.parallel import paropen
from gpaw.atom.basis import BasisMaker
from gpaw import GPAW, Mixer
from gpaw.response.grid_chi import CHI
from gpaw.mpi import serial_comm, rank, size
from gpaw.utilities import devnull


if rank != 0:
    sys.stdout = devnull 

GS = 1
GS2 = 1
OpticalLimit = 1

nband = 60


if GS:

    # Gs calculation one
    basis = BasisMaker('C').generate(2, 1) # dzp

    kpts = (20,20,7)
    a=1.42
    c=3.355

    # AB stack
    atoms = Atoms('C4',[
                  (1/3.0,1/3.0,0),
                  (2/3.0,2/3.0,0),
                  (0.   ,0.   ,0.5),
                  (1/3.0,1/3.0,0.5)
                  ],
                  pbc=(1,1,1))
    atoms.set_cell([(sqrt(3)*a/2.0,3/2.0*a,0),
                    (-sqrt(3)*a/2.0,3/2.0*a,0),
                    (0.,0.,2*c)],
                   scale_atoms=True)
    
    calc = GPAW(xc='LDA',
                kpts=kpts,
                h=0.2,
                basis='dzp',
                nbands=nband+10,
                convergence={'bands':nband},
                eigensolver='cg',
                symmetry='off',
                mixer=Mixer(0.05, 3, weight=100.0),
                width=0.05, txt='out.txt')
    
    atoms.set_calculator(calc)
#    view(atoms)

    atoms.get_potential_energy()
    calc.write('graphite1.gpw','all')


if GS2:

    # Gs calculation shifted kpoint
    calc = GPAW('graphite1.gpw')
    q = np.array([0.00001, 0., 0.])
    kpts = q + calc.get_ibz_k_points()
    assert calc.get_ibz_k_points().shape[0] == calc.get_bz_k_points().shape[0]

    atoms = calc.atoms
    calc = GPAW(xc='LDA',
                kpts=kpts,
                h = 0.2,
                basis='dzp',
                symmetry='off',
                nbands=nband+10,
                convergence={'bands':nband},
                eigensolver = 'cg',
                mixer=Mixer(0.05, 3, weight=100.0),
                width=0.05, txt = 'out2.txt')
    
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    calc.write('graphite2.gpw','all')    

    
if OpticalLimit:
    
    calc1 = GPAW('graphite1.gpw',communicator=serial_comm)
    calc2 = GPAW('graphite2.gpw',communicator=serial_comm)
                
    dw = 0.1  
    wmax = 40.

    q = np.array([0.00001, 0., 0.])
    chi = CHI()
    chi.nband = nband
    chi.initialize((calc1,calc2), q, wmax, dw, eta=0.2, Ecut = 40)
    chi.periodic()
    
    chi.get_absorption_spectrum('graphite_absorption')
    chi.check_sum_rule()

