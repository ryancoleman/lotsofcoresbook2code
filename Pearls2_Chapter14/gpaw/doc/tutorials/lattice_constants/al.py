import numpy as np
from ase.lattice import bulk
from gpaw import GPAW, PW, FermiDirac, MethfesselPaxton

a0 = 4.04
al = bulk('Al', 'fcc', a=a0)
cell0 = al.cell

for ecut in range(200, 501, 50):
    al.calc = GPAW(mode=PW(ecut),
                   xc='PBE',
                   kpts=(8, 8, 8),
                   parallel={'band': 1},
                   basis='dzp',
                   txt='Al-%d.txt' % ecut)
    for eps in np.linspace(-0.02, 0.02, 5):
        al.cell = (1 + eps) * cell0
        al.get_potential_energy()

al.calc.set(mode=PW(400))
for k in range(4, 17):
    al.calc.set(kpts=(k, k, k),
                txt='Al-%02d.txt' % k)
    for eps in np.linspace(-0.02, 0.02, 5):
        al.cell = (1 + eps) * cell0
        al.get_potential_energy()
