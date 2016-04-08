import os
from math import pi, cos, sin
from ase import Atom, Atoms
from ase.parallel import rank, barrier
from gpaw import GPAW, FermiDirac
from gpaw.test import equal, gen
from gpaw.xas import  * #XAS, RecursionMethod
import numpy as np
from gpaw import setup_paths
#setup_paths.insert(0, '.')


# Generate setup for oxygen with half a core-hole:
gen('Si', name='hch1s', corehole=(1, 0, 0.5))

a = 5.43095
b = a / 2
c = b / 2
d = b + c

si_nonortho = Atoms([Atom('Si', (0, 0, 0)),
            Atom('Si', (a/4, a/4, a/4))],
            cell=[(a/2,a/2,0),(a/2,0,a/2),(0,a/2,a/2)],
            pbc=True)

# calculation with full symmetry
calc = GPAW(nbands=-10,
            h=0.25,
            kpts=(2,2,2),
            occupations=FermiDirac(width=0.05),
            setups={0: 'hch1s'})

si_nonortho.set_calculator(calc)
e = si_nonortho.get_potential_energy()
niter = calc.get_number_of_iterations()
calc.write('si_nonortho_xas_sym.gpw')


# calculation without any symmetry
calc = GPAW(nbands=-10,
            h=0.25,
            kpts=(2,2,2),
            occupations=FermiDirac(width=0.05),
            setups={0: 'hch1s'},
            symmetry='off')

si_nonortho.set_calculator(calc)
e = si_nonortho.get_potential_energy()
niter = calc.get_number_of_iterations()
calc.write('si_nonortho_xas_nosym.gpw')


# restart from file
calc1 = GPAW('si_nonortho_xas_sym.gpw')
calc2 = GPAW('si_nonortho_xas_nosym.gpw')
if mpi.size == 1:
    xas1 = XAS(calc1)
    x, y1 = xas1.get_spectra()
    xas2 = XAS(calc2)
    x2, y2 = xas2.get_spectra(E_in=x)

    assert (np.sum(abs(y1-y2)[0,:500] ** 2 ) < 5e-9)
    assert (np.sum(abs(y1-y2)[1,:500] ** 2 ) < 5e-9)
    assert (np.sum(abs(y1-y2)[2,:500] ** 2 ) < 5e-9)

#else:
#    x = np.linspace(0, 10, 50)
#
#r1 = RecursionMethod(calc1)
#r1.run(40)
#
#r2 = RecursionMethod(calc1)
#r2.run(40)
#
#if mpi.size == 1:                                                                                                                                                                                                  
#    y1 = r1.get_spectra(x)
#    y2 = r2.get_spectra(x)
#
#    assert (np.sum(abs(y1-y2)[0,:] ** 2 ) < 5e-9)
#    assert (np.sum(abs(y1-y2)[1,:] ** 2 ) < 5e-9)
#    assert (np.sum(abs(y1-y2)[2,:] ** 2 ) < 5e-9)

