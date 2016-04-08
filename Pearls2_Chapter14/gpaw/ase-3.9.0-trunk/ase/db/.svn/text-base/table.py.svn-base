from __future__ import print_function

import numpy as np

from ase.data import atomic_masses
from ase.db.core import float_to_time_string, now, dict2constraint
from ase.utils import hill


all_columns = ['id', 'age', 'user', 'formula', 'calculator',
               'energy', 'fmax', 'pbc', 'volume',
               'charge', 'mass', 'smax', 'magmom']


def plural(n, word):
    if n == 1:
        return '1 ' + word
    return '%d %ss' % (n, word)

    
def cut(txt, length):
    if len(txt) <= length or length == 0:
        return txt
    return txt[:length - 3] + '...'


def cutlist(lst, length):
    if len(lst) <= length or length == 0:
        return lst
    return lst[:9] + ['... ({0} more)'.format(len(lst) - 9)]

    
def dict2forces(d):
    forces = d.get('forces')
    if forces is None:
        return None
        
    constraints = [dict2constraint(c) for c in d.get('constraints', [])]
    if constraints:
        forces = forces.copy()
        for constraint in constraints:
            constraint.adjust_forces(d.positions, forces)
            
    return forces

    
class Table:
    def __init__(self, connection, verbosity=1, cut=35):
        self.connection = connection
        self.verbosity = verbosity
        self.cut = cut
        self.rows = []
        self.columns = None
        self.id = None
        self.right = None
        self.keys = None
        
    def select(self, query, columns, sort, limit, offset):
        self.limit = limit
        self.offset = offset
        
        if sort != 'id':
            limit = 0
            offset = 0
           
        self.rows = [Row(d, columns)
                     for d in self.connection.select(
                         query, verbosity=self.verbosity,
                         limit=limit, offset=offset)]

        delete = set(range(len(columns)))
        for row in self.rows:
            for n in delete.copy():
                if row.values[n] is not None:
                    delete.remove(n)
        delete = sorted(delete, reverse=True)
        for row in self.rows:
            for n in delete:
                del row.values[n]
                
        self.columns = list(columns)
        for n in delete:
            del self.columns[n]
            
        if sort != 'id':
            reverse = sort[0] == '-'
            n = self.columns.index(sort.lstrip('-'))
            
            def key(row):
                x = row.values[n]
                return (x is None, x)
                
            self.rows = sorted(self.rows, key=key, reverse=reverse)
            
            if self.limit:
                self.rows = self.rows[self.offset:self.offset + self.limit]
                
    def format(self, subscript=None):
        right = set()
        allkeys = set()
        for row in self.rows:
            numbers = row.format(self.columns, subscript)
            right.update(numbers)
            allkeys.update(row.dct.key_value_pairs)
            
        right.add('age')
        self.right = [column in right for column in self.columns]
        
        self.keys = sorted(allkeys)

    def write(self):
        self.format()
        L = [[len(s) for s in row.strings]
             for row in self.rows]
        L.append([len(c) for c in self.columns])
        N = np.max(L, axis=0)

        fmt = '{0:{align}{width}}'
        print('|'.join(fmt.format(c, align='<>'[a], width=w)
                       for c, a, w in zip(self.columns, self.right, N)))
        for row in self.rows:
            print('|'.join(fmt.format(c, align='<>'[a], width=w)
                           for c, a, w in
                           zip(row.strings, self.right, N)))

        if self.verbosity == 0:
            return
            
        print('Rows:', len(self.rows), end='')
        if self.limit and len(self.rows) == self.limit:
            print(' (limited to first {0})'.format(self.limit))
        else:
            print()

        if self.keys:
            print('Keys:', ', '.join(cutlist(self.keys, self.cut)))
            
    def write_csv(self):
        print(', '.join(self.columns))
        for row in self.rows:
            print(', '.join(str(val) for val in row.values))

            
class Row:
    def __init__(self, dct, columns):
        self.dct = dct
        self.values = None
        self.strings = None
        self.more = False
        self.set_columns(columns)
        if 'key_value_pairs' not in dct:
            dct['key_value_pairs'] = {}
        
    def set_columns(self, columns):
        self.values = []
        for c in columns:
            f = getattr(self, c, None)
            if f is None:
                value = getattr(self.dct, c, None)
            else:
                try:
                    value = f(self.dct)
                except (AttributeError, TypeError):
                    value = None
            self.values.append(value)
            
    def toggle(self):
        self.more = not self.more
        
    def format(self, columns, subscript=None):
        self.strings = []
        numbers = set()
        for value, column in zip(self.values, columns):
            if column == 'formula' and subscript:
                value = subscript.sub(r'<sub>\1</sub>', value)
            elif isinstance(value, int):
                value = str(value)
                numbers.add(column)
            elif isinstance(value, float):
                numbers.add(column)
                value = '{0:.3f}'.format(value)
            elif value is None:
                value = ''
            self.strings.append(value)
        
        return numbers
        
    def age(self, d):
        return float_to_time_string(now() - d.ctime)

    def formula(self, d):
        return hill(d.numbers)

    def volume(self, d):
        return abs(np.linalg.det(d.cell))

    def pbc(self, d):
        return ''.join('-P'[p] for p in d.pbc)

    def fmax(self, d):
        forces = dict2forces(d)
        return (forces**2).sum(1).max()**0.5

    def mass(self, d):
        if 'masses' in d:
            return d.masses.sum()
        return atomic_masses[d.numbers].sum()

    def smax(self, d):
        return (d.stress**2).max()**0.5
