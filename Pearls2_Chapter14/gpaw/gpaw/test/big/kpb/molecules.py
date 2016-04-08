"""Test exact exchange for 20 small molecules.

Compare results to::

  S. Kurth, J. P. Perdew, and P. Blaha
  Molecular and Soild-State Tests of Density Functional
  Approximations: LSD, GGAs, and Meta-GGAs
  International Journal of Quantum Chemistry, Vol. 75, 889-909, 1999

"""

import ase.db
from ase import Atoms
from ase.structure import molecule
from ase.data.g2_1_ref import diatomic, ex_atomization

from gpaw import GPAW, PW
from gpaw.xc.exx import EXX


# Experimental bondlengths:
bondlengths = {'H2': 0.741,
               'OH': 0.970,
               'HF': 0.9168,
               'NO': 1.154,
               'P2': 1.893}
bondlengths.update((name, d[0]) for name, d in diatomic.items())

extra = {
   'CH4': ('CH4', [(0.0000, 0.0000, 0.0000),
                   (0.6276, 0.6276, 0.6276),
                   (0.6276, -0.6276, -0.6276),
                   (-0.6276, 0.6276, -0.6276),
                   (-0.6276, -0.6276, 0.6276)]),
   'NH3': ('NH3', [(0.0000, 0.0000, 0.0000),
                   (0.0000, -0.9377, -0.3816),
                   (0.8121, 0.4689, -0.3816),
                   (-0.8121, 0.4689, -0.3816)]),
   'H2O': ('OH2', [(0.0000, 0.0000, 0.1173),
                   (0.0000, 0.7572, -0.4692),
                   (0.0000, -0.7572, -0.4692)]),
   'C2H2': ('C2H2', [(0.0000, 0.0000, 0.6013),
                     (0.0000, 0.0000, -0.6013),
                     (0.0000, 0.0000, 1.6644),
                     (0.0000, 0.0000, -1.6644)]),
   'C2H4': ('C2H4', [(0.0000, 0.0000, 0.6695),
                     (0.0000, 0.0000, -0.6695),
                     (0.0000, 0.9289, 1.2321),
                     (0.0000, -0.9289, 1.2321),
                     (0.0000, 0.9289, -1.2321),
                     (0.0000, -0.9289, -1.2321)]),
   'HCN': ('CHN', [(0.0000, 0.0000, 0.0000),
                   (0.0000, 0.0000, 1.0640),
                   (0.0000, 0.0000, -1.1560)]),
   'Be2': ('Be2', [(0.0000, 0.0000, 0.0000),
                   (0.0000, 0.0000, 2.460)])}


c = ase.db.connect('results.db')

for name in ex_atomization.keys() + 'H Li Be B C N O F Cl P'.split():
    id = c.reserve(name=name)
    if id is None:
        continue
        
    if name in extra:
        a = Atoms(*extra[name])
    else:
        a = molecule(name)
        if name in bondlengths:
            a.set_distance(0, 1, bondlengths[name])
    a.cell = [11, 12, 13]
    a.center()
   
    a.calc = GPAW(xc='PBE', mode=PW(500), txt=name + '.txt', dtype=complex)
    a.get_potential_energy()
    
    exx = EXX(a.calc)
    exx.calculate()
    eexx = exx.get_total_energy()

    c.write(a, name=name, exx=eexx)
    del c[id]
    