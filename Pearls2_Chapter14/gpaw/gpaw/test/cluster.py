from __future__ import print_function
import numpy as np

from ase import Atoms, Atom
from ase.parallel import barrier, rank, size
from gpaw.cluster import Cluster
from gpaw.test import equal
from ase.structure import molecule
from math import pi, sqrt

R = 2.0
CO = Atoms([Atom('C', (1, 0, 0)), Atom('O', (1, 0, R))])

CO.rotate('y', pi/2)
equal(CO.positions[1, 0], R, 1e-10)

# translate
CO.translate(-CO.get_center_of_mass())
p = CO.positions.copy()
for i in range(2):
    equal(p[i, 1], 0, 1e-10)
    equal(p[i, 2], 0, 1e-10)

# rotate the nuclear axis to the direction (1,1,1)
CO.rotate(p[1] - p[0], (1, 1, 1))
q = CO.positions.copy()
for c in range(3):
    equal(q[0, c], p[0, 0] / sqrt(3), 1e-10)
    equal(q[1, c], p[1, 0] / sqrt(3), 1e-10)

# minimal box
b=4.0
CO = Cluster([Atom('C', (1, 0, 0)), Atom('O', (1, 0, R))])
CO.minimal_box(b)
cc = CO.get_cell() 
for c in range(3):
    width = 2*b
    if c==2:
        width += R
    equal(cc[c, c], width, 1e-10)

# minimal box, ensure multiple of 4
h = .13
b = [2, 3, 4]
CO.minimal_box(b, h=h)
cc = CO.get_cell() 
for c in range(3):
#    print "cc[c,c], cc[c,c] / h % 4 =", cc[c, c], cc[c, c] / h % 4
    for a in CO:
        print(a.symbol, b[c], a.position[c], cc[c, c] - a.position[c])
        assert(a.position[c] > b[c])
    equal(cc[c, c] / h % 4, 0.0, 1e-10)

# .............................................
# connected atoms
assert(len(CO.find_connected(0, 1.1 * R)) == 2)
assert(len(CO.find_connected(0, 0.9 * R)) == 1)

H2O = Cluster(molecule('H2O'))
assert (len(H2O.find_connected(0)) == 3)
assert (len(H2O.find_connected(0, scale=0.9)) == 1)

# .............................................
# I/O
fxyz='CO.xyz'
fpdb='CO.pdb'

cell = [2.,3.,R+2.]
CO.set_cell(cell, scale_atoms=True)
barrier()
CO.write(fxyz)
barrier()
CO_b = Cluster(filename=fxyz)
assert(len(CO) == len(CO_b))
#for a, b in zip(cell, CO_b.get_cell().diagonal()):
#    assert(a == b)
offdiagonal = CO_b.get_cell().sum() - CO_b.get_cell().diagonal().sum()
assert(offdiagonal == 0.0)
 
barrier()
CO.write(fpdb)

# read xyz files with additional info
read_with_additional = True
if read_with_additional:
    if rank == 0:
        f = open(fxyz, 'w')
        print("""2

C 0 0 0. 1 2 3
O 0 0 1. 6. 7. 8.""", file=f)
        f.close()

    barrier()

    CO = Cluster(filename=fxyz)
