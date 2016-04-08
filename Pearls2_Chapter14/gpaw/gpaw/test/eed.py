import os

from ase import Atom 
from gpaw import GPAW
from gpaw.cluster import Cluster
from gpaw.test import equal
from gpaw.analyse.eed import ExteriorElectronDensity

fwfname='H2_kpt441_wf.gpw'
txt = None
#txt = '-'

# write first if needed
try:
    c = GPAW(fwfname, txt=txt)
    s = c.get_atoms()
except:
    s = Cluster([Atom('H'), Atom('H', [0,0,1])], pbc=[1,1,0])
    s.minimal_box(3.)
    c = GPAW(xc='PBE', h=.3, kpts=(4,4,1),
             convergence={'density':1e-3, 'eigenstates':1e-3})
    c.calculate(s)
    c.write(fwfname, 'all')

eed = ExteriorElectronDensity(c.wfs.gd, s)
eed.write_mies_weights(c.wfs)
