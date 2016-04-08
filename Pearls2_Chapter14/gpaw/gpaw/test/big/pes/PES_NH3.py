# Takes 30 mins with 8 cores in domain-parallel (memory ~500mb, max 1.5gb)

import numpy as np
from ase import Atom, Atoms
from gpaw import GPAW, MixerDif, mpi
from gpaw.lrtddft import LrTDDFT
from gpaw.pes.dos import DOSPES
from gpaw.pes.tddft import TDDFTPES

atoms = Atoms([Atom('N', (5.915080000001099, 5.584099935483791, 5.517093011559592)),
               Atom('H', (5.915080000000824, 6.534628192918412, 5.154213005717621)),
               Atom('H', (6.739127072810322, 5.112845817889995, 5.151626633318725)),
               Atom('H', (5.091032927190902, 5.112845817890320, 5.151626633319492))])
atoms.center(vacuum=8)

h = 0.15
N_c = np.round(atoms.get_cell().diagonal() / h / 16) * 16

m_c = GPAW(gpts=N_c, nbands=4, mixer=MixerDif(0.1, 5, weight=100.0),
           #convergence={'bands':10},
           parallel={'domain': mpi.size},
           xc='PBE', txt='NH3-m.txt', spinpol=True)

m = atoms.copy()
m.set_initial_magnetic_moments([-1,1,1,-1])
m.set_calculator(m_c)
m.get_potential_energy()

d_c = GPAW(gpts=N_c, nbands=16, mixer=MixerDif(0.1, 5, weight=100.0),
           convergence={'bands':10},
           parallel={'domain': mpi.size},
           xc='PBE', txt='NH3-d.txt', spinpol=True)

d = atoms.copy()
d.set_initial_magnetic_moments([-1, 0.5, 0.5, -0.5])
d_c.set(charge=1)
d.set_calculator(d_c)
d.get_potential_energy()

istart=0 # band index of the first occ. band to consider
jend=15  # band index of the last unocc. band to consider
d_lr = LrTDDFT(d_c, xc='PBE', nspins=2 , istart=istart, jend=jend)
d_lr.diagonalize()

pes = TDDFTPES(m_c, d_lr, d_c)
pes.save_folded_pes('NH3-td.dat', folding=None)

pes = DOSPES(m_c, d_c, shift=True)
pes.save_folded_pes('NH3-dos.dat', folding=None)

