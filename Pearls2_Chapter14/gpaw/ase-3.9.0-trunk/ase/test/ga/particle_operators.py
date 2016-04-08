from ase.cluster import Icosahedron
from ase.ga.particle_crossovers import CutSpliceCrossover
from random import shuffle

ico1 = Icosahedron('Cu', 5)
ico1.info['confid'] = 1
ico2 = Icosahedron('Ni', 5)
ico2.info['confid'] = 2

op = CutSpliceCrossover({(28, 29): 2.0, (28, 28): 2.0, (29, 29): 2.0},
                        keep_composition=False)
a3, desc = op.get_new_individual([ico1, ico2])

assert len(set(a3.get_chemical_symbols())) == 2

ico1.numbers[:147] = [28] * 147
shuffle(ico1.numbers)
ico2.numbers[:162] = [29] * 162
shuffle(ico2.numbers)
op = CutSpliceCrossover({(28, 29): 2.0, (28, 28): 2.0, (29, 29): 2.0})
a3, desc = op.get_new_individual([ico1, ico2])

# assert a3.get_chemical_formula() == 'Cu162Ni147'

from ase.ga.particle_mutations import COM2surfPermutation
# from ase.ga.particle_mutations import RandomPermutation
# from ase.ga.particle_mutations import Poor2richPermutation
# from ase.ga.particle_mutations import Rich2poorPermutation

op = COM2surfPermutation()
a3, desc = op.get_new_individual([ico1])
a3.info['confid'] = 3

assert a3.get_chemical_formula() == 'Cu162Ni147'

atomic_conf = op.get_atomic_configuration(a3, elements=['Cu'])[-4:]
cu3 = len([item for sublist in atomic_conf for item in sublist])
a4, desc = op.get_new_individual([a3])
atomic_conf = op.get_atomic_configuration(a4, elements=['Cu'])[-4:]
cu4 = len([item for sublist in atomic_conf for item in sublist])

assert abs(cu4 - cu3) == 1
