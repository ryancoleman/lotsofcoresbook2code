from __future__ import print_function
from gpaw import GPAW

print('state    LDA        PBE')

for name in ['ferro', 'anti', 'non']:

    calc = GPAW(name + '.gpw', txt=None)
    atoms = calc.get_atoms()
    eLDA = atoms.get_potential_energy()
    deltaxc = calc.get_xc_difference('PBE')
    ePBE = eLDA + deltaxc

    if name == 'ferro':
        eLDA0 = eLDA
        ePBE0 = ePBE

    eLDA -= eLDA0
    ePBE -= ePBE0
    print(('%-5s: %7.3f eV %7.3f eV' % (name, eLDA, ePBE)))
