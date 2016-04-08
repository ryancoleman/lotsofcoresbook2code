from __future__ import print_function
from ase.structure import molecule
from gpaw import GPAW

a = 8.0
h = 0.2

energies = {}
resultfile = open('results-%.2f.txt' % h, 'w')

for name in ['H2O', 'H', 'O']:
    system = molecule(name)
    system.set_cell((a, a, a))
    system.center()

    calc = GPAW(h=h,
                txt='gpaw-%s-%.2f.txt' % (name, h))
    if name == 'H' or name == 'O':
        calc.set(hund=True)

    system.set_calculator(calc)

    energy = system.get_potential_energy()
    energies[name] = energy
    print(name, energy, file=resultfile)

e_atomization = energies['H2O'] - 2 * energies['H'] - energies['O']
print(e_atomization, file=resultfile)
