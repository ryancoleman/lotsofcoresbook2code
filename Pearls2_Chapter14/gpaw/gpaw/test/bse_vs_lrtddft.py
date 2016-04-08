from __future__ import print_function
import numpy as np
from ase import Atom, Atoms
from ase.units import Hartree
from gpaw.mpi import size
from gpaw import GPAW
from gpaw.response.bse import BSE

GS = 1
bse = 1
casida = 1
compare = 1

if GS:
    d = 2.89
    cluster = Atoms([Atom('Na', (0, 0, 0)),
                     Atom('Na', (0, 0, d)),
                     ], pbc=True)
    cluster.set_cell((15.,15.,18.), scale_atoms=False)
    cluster.center()
    calc = GPAW(h=0.3, nbands=8, setups={'Na': '1'})

    cluster.set_calculator(calc)
    cluster.get_potential_energy()
    calc.write('Na2.gpw','all')

if bse:    
    bse = BSE('Na2.gpw',
              w=np.linspace(0,15,151),
              nv=[0,8],
              nc=[0,8],
              mode='RPA',
              coupling=True,
              q=np.array([0,0,0.0001]),
              optical_limit=True,
              ecut=50.,
              nbands=8)
    bse.initialize()
    H_SS = bse.calculate()
    bse.diagonalize(H_SS)
    
    w = np.real(bse.w_S) * Hartree
    print(np.shape(w))
    energies = np.sort(w)[len(w)/2:]
    print('BSE:', energies)

if casida:
    from gpaw.lrtddft import LrTDDFT
    from gpaw.lrtddft import photoabsorption_spectrum

    calc = GPAW('Na2.gpw',txt=None)

    lr = LrTDDFT(calc, xc=None, istart=0, jend=7, nspins=1) 
    lr.diagonalize()
    photoabsorption_spectrum(lr, 'Na2_spectrum.dat', width=0.05)   

    energies_lrtddft = lr.get_energies() * Hartree
    print('lrTDDFT:', energies_lrtddft)
    
if compare:
    assert (np.abs(energies - energies_lrtddft)).max() < 3*1e-3
