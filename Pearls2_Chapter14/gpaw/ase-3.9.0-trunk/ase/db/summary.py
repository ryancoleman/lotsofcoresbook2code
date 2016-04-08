from __future__ import print_function

from ase.db.table import dict2forces
from ase.data import atomic_masses, chemical_symbols
from ase.db.core import float_to_time_string, now
from ase.utils import hill

import numpy as np


class Summary:
    def __init__(self, dct, subscript=None):
        self.dct = dct
        
        self.cell = [['{0:.3f}'.format(a) for a in axis] for axis in dct.cell]
        
        forces = dict2forces(dct)
        if forces is None:
            fmax = None
            self.forces = None
        else:
            fmax = (forces**2).sum(1).max()**0.5
            N = len(forces)
            self.forces = []
            for n, f in enumerate(forces):
                if n < 5 or n >= N - 5:
                    f = tuple('{0:10.3f}'.format(x) for x in f)
                    symbol = chemical_symbols[dct.numbers[n]]
                    self.forces.append((n, symbol) + f)
                elif n == 5:
                    self.forces.append((' ...', '',
                                        '       ...',
                                        '       ...',
                                        '       ...'))
                    
        self.stress = dct.get('stress')
        if self.stress is not None:
            self.stress = ', '.join('{0:.3f}'.format(s) for s in self.stress)
            
        if 'masses' in dct:
            mass = dct.masses.sum()
        else:
            mass = atomic_masses[dct.numbers].sum()
            
        formula = hill(dct.numbers)
        if subscript:
            formula = subscript.sub(r'<sub>\1</sub>', formula)
            
        table = [
            ('id', '', dct.id),
            ('age', '', float_to_time_string(now() - dct.ctime, True)),
            ('formula', '', formula),
            ('user', '', dct.user),
            ('calculator', '', dct.get('calculator')),
            ('energy', 'eV', dct.get('energy')),
            ('fmax', 'eV/Ang', fmax),
            ('charge', '|e|', dct.get('charge')),
            ('mass', 'au', mass),
            ('unique id', '', dct.unique_id),
            ('volume', 'Ang^3', abs(np.linalg.det(dct.cell)))]
        self.table = [(name, unit, value) for name, unit, value in table
                      if value is not None]

        if 'key_value_pairs' in dct:
            self.key_value_pairs = sorted(dct.key_value_pairs.items())
        else:
            self.key_value_pairs = None

        self.dipole = dct.get('dipole')
        if self.dipole is not None:
            self.dipole = ', '.join('{0:.3f}'.format(d) for d in self.dipole)
        
        self.data = dct.get('data')
        if self.data:
            self.data = ', '.join(self.data.keys())
            
        self.constraints = dct.get('constraints')
        if self.constraints:
            self.constraints = ', '.join(d['name'] for d in self.constraints)
        
    def write(self):
        dct = self.dct
        
        width = max(len(name) for name, unit, value in self.table)
        print('{0:{width}}|unit  |value'.format('name', width=width))
        for name, unit, value in self.table:
            print('{0:{width}}|{1:6}|{2}'.format(name, unit, value,
                                                 width=width))

        print('\nUnit cell in Ang:')
        print('axis|periodic|          x|          y|          z')
        c = 1
        for p, axis in zip(dct.pbc, self.cell):
            print('   {0}|     {1}|{2[0]:>11}|{2[1]:>11}|{2[2]:>11}'.format(
                c, [' no', 'yes'][p], axis))
            c += 1
            
        if self.key_value_pairs:
            print('\nKey-value pairs:')
            width = max(len(key) for key, value in self.key_value_pairs)
            for key, value in self.key_value_pairs:
                print('{0:{width}}|{1}'.format(key, value, width=width))
                
        if self.forces:
            print('\nForces in ev/Ang:')
            for f in self.forces:
                print('{0:4}|{1:2}|{2}|{3}|{4}'.format(*f))

        if self.stress:
            print('\nStress tensor (xx, yy, zz, zy, zx, yx) in eV/Ang^3:')
            print('   ', self.stress)

        if self.dipole:
            print('\nDipole moment in e*Ang: ({0})'.format(self.dipole))
        
        if self.constraints:
            print('\nConstraints:', self.constraints)
            
        if self.data:
            print('\nData:', self.data)
