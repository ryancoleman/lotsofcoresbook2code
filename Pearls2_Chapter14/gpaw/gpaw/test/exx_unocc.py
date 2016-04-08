from __future__ import print_function

from ase import Atoms
from ase.parallel import parprint
from ase.utils.timing import Timer

from gpaw import GPAW
from gpaw.test import equal
from gpaw.xc.hybrid import HybridXC

timer = Timer()

loa = Atoms('Be2',
            [(0, 0, 0), (2.45, 0, 0)],
            cell=[5.9, 4.8, 5.0])
loa.center()

txt = None
xc = 'PBE0'
nbands = 4

unocc = True
load = False

# usual calculation
fname = 'Be2.gpw'
if not load:
    xco = HybridXC(xc)
    cocc = GPAW(h=0.3,
                eigensolver='rmm-diis',
                xc=xco,
                nbands=nbands,
                convergence={'eigenstates': 1e-4},
                txt=txt)
    cocc.calculate(loa)
else:
    cocc = GPAW(fname)
    cocc.converge_wave_functions()
fo_n = 1. * cocc.get_occupation_numbers()
eo_n = 1. * cocc.get_eigenvalues()

if unocc:
    # apply Fock opeartor also to unoccupied orbitals
    xcu = HybridXC(xc, unocc=True)
    cunocc = GPAW(h=0.3,
                  eigensolver='rmm-diis',
                  xc=xcu,
                  nbands=nbands,
                  convergence={'eigenstates': 1e-4},
                  txt=txt)
    cunocc.calculate(loa)

    parprint('     HF occ          HF unocc      diff')
    parprint('Energy %10.4f   %10.4f %10.4f' %
             (cocc.get_potential_energy(),
              cunocc.get_potential_energy(),
              cocc.get_potential_energy() - cunocc.get_potential_energy()
              ))
    equal(cocc.get_potential_energy(),
          cunocc.get_potential_energy(), 1.e-4)

    fu_n = cunocc.get_occupation_numbers()
    eu_n = cunocc.get_eigenvalues()

    parprint('Eigenvalues:')
    for eo, fo, eu, fu in zip(eo_n, fo_n, eu_n, fu_n):
        parprint('%8.4f %5.2f   %8.4f %5.2f  %8.4f' %
                 (eo, fo, eu, fu, eu - eo))
        if fo > 0.01:
            equal(eo, eu, 3.5e-4)
