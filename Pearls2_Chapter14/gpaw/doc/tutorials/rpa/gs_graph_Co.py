from __future__ import print_function
import numpy as np
from ase.visualize import view
from ase.dft.kpoints import monkhorst_pack
from ase.parallel import paropen
from ase.lattice.surface import *
from gpaw import GPAW, PW, FermiDirac, MixerSum
from gpaw.xc.exx import EXX

kpts = monkhorst_pack((16, 16, 1))
kpts += np.array([1/32., 1/32., 0])

a = 2.51 # Lattice parameter of Co
slab = hcp0001('Co', a=a, c=4.07, size=(1,1,4))
pos = slab.get_positions()
cell = slab.get_cell()
cell[2,2] = 20. + pos[-1, 2]
slab.set_cell(cell)
slab.set_initial_magnetic_moments([0.7, 0.7, 0.7, 0.7])

ds = [1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 5.0, 6.0, 10.0]

for d in ds:
    pos = slab.get_positions()
    add_adsorbate(slab, 'C', d, position=(pos[3,0], pos[3,1]))
    add_adsorbate(slab, 'C', d, position=(cell[0,0]/3 + cell[1,0]/3,
                                          cell[0,1]/3 + cell[1,1]/3))
    #view(slab)
    calc = GPAW(xc='PBE',
                eigensolver='cg',
                mode=PW(600),
                kpts=kpts,
                occupations=FermiDirac(width=0.01),
                mixer=MixerSum(beta=0.1, nmaxold=5, weight=50.0),
                convergence={'density': 1.e-6},
                maxiter=300,
                parallel={'domain': 1,
                          'band': 1},
                txt='gs_%s.txt' % d)
    slab.set_calculator(calc)
    E = slab.get_potential_energy()
    exx = EXX(calc, txt='exx_%s.txt' % d)
    exx.calculate()
    E_hf = exx.get_total_energy()

    calc.diagonalize_full_hamiltonian()
    calc.write('gs_%s.gpw' % d, mode='all')

    f = paropen('hf_acdf.dat', 'a')
    print(d, E, E_hf, file=f)
    f.close()
    
    del slab[-2:]
