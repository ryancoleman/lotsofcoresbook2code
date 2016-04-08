from __future__ import print_function
from gpaw import restart

# J.Phys.: Condens. Matter 18 (2006) 41-54
pw91vasp = {'NO':     -0.95,
            'N2':      0.00,
            'O2':      0.00,
            'NRu001': -0.94,
            'ORu001': -2.67,
            'Ru001':   0.00}

pbe = {}
pw91 = {}
for name in ['NO', 'O2', 'N2', 'Ru001', 'NRu001', 'ORu001']:
    a, calc = restart(name, txt=None)
    pbe[name] = a.get_potential_energy()
    pw91[name] = pbe[name] + calc.get_xc_difference('PW91')

for data, text in [(pbe, 'PBE'),
                   (pw91, 'PW91 (non-selfconsitent)'),
                   (pw91vasp, 'PW91 (VASP)')]:
    print(('%22s %.3f %.3f %.3f' %
           (text,
            data['NRu001'] - data['Ru001'] - data['N2'] / 2,
            data['ORu001'] - data['Ru001'] - data['O2'] / 2,
            data['NO'] - data['N2'] / 2 - data['O2'] / 2)))

