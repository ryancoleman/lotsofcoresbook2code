from __future__ import print_function

import sys
from math import sqrt, pi

import numpy as np

from ase import Atoms
from ase.units import Bohr, Hartree
from ase.parallel import rank
from gpaw import GPAW

from gpaw.external_potential import ConstantElectricField

###

field   = 0.01
dx      = 2.0
vac     = 3.0
nsteps  = 3

nbands  = 4

a0      = dx+2*vac
b0      = 2*vac

debug   = False

if debug:
    txt = 'gpaw.out'
else:
    txt = None

###

a  = Atoms('Be', positions=[ [ b0/2, b0/2, a0/2  ] ], cell=[ b0, b0, a0 ])

z_list = np.linspace(a0/2-dx/2, a0/2+dx/2, nsteps)
    
if True:
    c = GPAW(
        h         = 0.30,
        width     = 0.0,
        nbands    = nbands,
        spinpol   = False,
        xc        = 'LDA',
        txt       = txt,
        )
    a.set_calculator(c)

    ###

    e_no_field = [ ]
    for z in z_list:
        if rank == 0 and debug:
            print(z)
    
        a[0].z = z

        e_no_field += [ a.get_potential_energy() ]

    e_no_field = np.array(e_no_field)

###

if True:
    c = GPAW(
        h         = 0.30,
        width     = 0.0,
        nbands    = nbands,
        spinpol   = False,
        xc        = 'LDA',
        txt       = txt,
        external  = ConstantElectricField(field * Bohr/Hartree)
        )
    a.set_calculator(c)


    e_with_field = [ ]
    for z in z_list:
        if rank == 0 and debug:
            print(z)
    
        a[0].z = z

        e_with_field += [ a.get_potential_energy() ]

    e_with_field = np.array(e_with_field)

###

np.savetxt("e.out", np.transpose( [ z_list, e_with_field-e_no_field, e_no_field, e_with_field ] ))

c1, c2 = np.polyfit(z_list, e_with_field-e_no_field, 1)

# 4*field because the charge of the nuclei is not considered
# in ExternalPotential
err    = abs(c1-a[0].number*field)

if rank == 0 and debug:
    print(c1)
    print(err)

assert err < 0.001
