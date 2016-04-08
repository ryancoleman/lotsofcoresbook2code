from __future__ import print_function
import sys
import re

import numpy as np

from ase.visualize import view
from ase.structure import molecule
from ase.units import Hartree
from gpaw import GPAW
from gpaw.mpi import rank, world
from gpaw.test import equal
from gpaw.gauss import Gauss
from gpaw.lrtddft import LrTDDFT, photoabsorption_spectrum
from gpaw.lrtddft.kssingle import KSSingles
from cStringIO import StringIO

L = 10.0
txt=None
txt='-'

N2 = molecule('N2')
N2.set_cell([L, L, L])
#N2.set_pbc(True)
N2.center()

try:
    calc = GPAW('N2_wfs.gpw', 
                txt=txt, 
                parallel={'domain': world.size})
    calc.converge_wave_functions()
except:
    calc = GPAW(h=0.25,
                nbands=-5,
                spinpol=True,
                xc='PBE',
                txt=txt,
                eigensolver='cg',
                parallel={'domain': world.size})
    N2.set_calculator(calc)
    E0 = N2.get_potential_energy()
    calc.write('N2_wfs.gpw', 'all')

# selections
for obj in [KSSingles, LrTDDFT]:
    # selection using state numbers
    el = obj(calc, istart=3, jend=6, txt=txt)
    if hasattr(obj, 'diagonalize'):
        el.diagonalize()
#    print "*************** obj, len(obj)", obj.__name__, len(el)
    assert len(el) == 8
    # selection using an energy range
    el = obj(calc, energy_range=8, txt=txt)
    if hasattr(obj, 'diagonalize'):
        el.diagonalize()
#    print "*************** obj, len(obj)", obj.__name__, len(el)
    assert len(el) == 4
    el = obj(calc, energy_range=11.5, txt=txt)
#    print "*************** obj, len(obj)", obj.__name__, len(el)
    if hasattr(obj, 'diagonalize'):
        el.diagonalize()
    assert len(el) == 18
    if hasattr(obj, 'diagonalize'):
        el.diagonalize(energy_range=8)
        assert len(el) == 4

lr = LrTDDFT(calc, nspins=2)
lr.write('lrtddft3.dat.gz')
lr.diagonalize()

world.barrier()

# This is done to test if writing and reading again yields the same result
lr2 = LrTDDFT('lrtddft3.dat.gz')
lr2.diagonalize()

# Unfortunately not all of the lrtddft code is parallel
if rank == 0:
    Epeak = 19.5# The peak we want to investigate (this is alone)
    Elist = np.asarray([lrsingle.get_energy() * Hartree for lrsingle in lr])
    n = np.argmin(np.abs(Elist - Epeak)) # Index of the peak

    E = lr[n].get_energy() * Hartree
    osz = lr[n].get_oscillator_strength()
    print('Original object        :', E, osz[0])

    # Test the output of analyse
    origstdout = sys.stdout
    sys.stdout = sio = StringIO()
    lr.analyse(n)
    s = sio.getvalue() 
    sys.stdout = origstdout
    match = re.findall(r'%i: E=([0-9]*\.[0-9]*) eV, f=([0-9]*\.[0-9]*)*' % n, s)
    Eanalyse = float(match[0][0])
    oszanalyse = float(match[0][1])
    print('From analyse           :', Eanalyse, oszanalyse)
    equal(E, Eanalyse, 1e-3)            # Written precision in analyse
    equal(osz[0], oszanalyse, 1e-3)

    E2 = lr2[n].get_energy() * Hartree
    osz2 = lr2[n].get_oscillator_strength()
    print('Written and read object:', E2, osz2[0])
    
    # Compare values of original and written/read objects   
    equal(E, E2, 1e-4)
    for i in range(len(osz)):
        equal(osz[i], osz2[i], 1.7e-4)

    width = 0.05
    photoabsorption_spectrum(lr, 
                             spectrum_file = 'lrtddft3-spectrum.dat', 
                             width = width)
    # We need to be able to check the heights in the spectrum
    weight = Gauss(width).get(0)

    spectrum = np.loadtxt('lrtddft3-spectrum.dat', usecols = (0, 1))
    idx = (spectrum[:, 0] >= E - 0.1) & (spectrum[:, 0] <= E + 0.1)
    peak = np.argmax(spectrum[idx, 1]) + np.nonzero(idx)[0][0]
    Espec = spectrum[peak, 0]
    oszspec = spectrum[peak, 1] / weight

    print('Values from spectrum   :', Espec, oszspec)
    # Compare calculated values with values written to file
    equal(E, Espec, 1e-2)           # The spectrum has a low sampling
    equal(osz[0], oszspec, 1e-2)
