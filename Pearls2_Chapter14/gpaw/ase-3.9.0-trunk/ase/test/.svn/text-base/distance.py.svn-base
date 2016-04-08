import itertools
import numpy as np

from ase import Atoms, Atom, visualize
from ase.utils.distance import distance

# artificial structure
org = Atoms([
        Atom('C', [-1.75072, 0.62689, 0.00000]),
        Atom('O', [0.58357, 2.71652, 0.00000]),
        Atom('P', [-5.18268, 1.36522, 0.00000]),
        Atom('N', [-1.86663, -0.77867, 2.18917]),
        Atom('S', [-1.80586, 0.20783 , -2.79331]),
        ])
#visualize.view(org)

maxdist = 3.0e-13

# translate
for dx in range(3, 10, 2):
    new = org.copy()
    new.translate([dx / np.sqrt(2), -dx / np.sqrt(2), 0])
    dist = distance(org, new, True)
    dist2 = distance(org, new, False)
    print 'translation', dx, '-> distance', dist
    assert dist < maxdist
    assert dist == dist2

# rotate
for axis in ['x', '-y', 'z', np.array([1, 1, 1] / np.sqrt(3))]:
    for rot in [20, 200]:
        new = org.copy()
        new.translate(-new.get_center_of_mass())
        new.rotate(axis, np.pi * rot / 180)
        dist = distance(org, new, True)
        dist2 = distance(org, new, False)
        print 'rotation', axis, ', angle', rot, '-> distance', dist
        assert dist < maxdist
        assert dist == dist2
    
if 0:
    # reflect
    new = Atoms()
    cm = org.get_center_of_mass()
    for a in  org:  
        new.append(Atom(a.symbol, -(a.position - cm)))
    dist = distance(org, new)
    print 'reflected -> distance', dist

# permute
for i, a in enumerate(org):
    if i < 3:
        a.symbol = 'H'

if hasattr(itertools, 'permutations'):
    for indxs in itertools.permutations(range(3)):
        new = org.copy()
        for c in range(3):
            new[c].position = org[indxs[c]].position
        dist = distance(org, new)
        print 'permutation', indxs, '-> distance', dist 
        assert dist < maxdist
