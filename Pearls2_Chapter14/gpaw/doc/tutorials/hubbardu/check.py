from __future__ import print_function
from ase.io import read


fd = open('gaps.csv', 'w')
gaps = []
for name in ['no_u', 'normalized_u', 'not_normalized_u']:
    n = read(name + '.txt')
    gap = n.calc.get_eigenvalues(spin=1)[1] - n.calc.get_eigenvalues(spin=0)[1]
    gaps.append(gap)
    print('%s, %.3f' % (name.replace('_', ' ').replace('u', 'U'), gap), file=fd)

assert abs(gaps[1] - gaps[0] - 6.0) < 0.8
