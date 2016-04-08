"""SQLite3 backend.

Versions:

1) Added 3 more columns.
2) Changed "user" to "username".
3) Now adding keys to keyword table and added an "information" table containing
   a version number.
4) Got rid of keywords.
"""

from __future__ import absolute_import, print_function
import sqlite3

import numpy as np

from ase.db.core import Database, ops, now, lock, parallel, invop
from ase.db.jsondb import encode, decode


VERSION = 4

init_statements = [
    """CREATE TABLE systems (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    unique_id TEXT UNIQUE,
    ctime REAL,
    mtime REAL,
    username TEXT,
    numbers BLOB,
    positions BLOB,
    cell BLOB,
    pbc INTEGER,
    initial_magmoms BLOB,
    initial_charges BLOB,
    masses BLOB,
    tags BLOB,
    momenta BLOB,
    constraints TEXT,
    calculator TEXT,
    calculator_parameters TEXT,
    energy REAL,
    free_energy REAL,
    forces BLOB,
    stress BLOB,
    dipole BLOB,
    magmoms BLOB,
    magmom BLOB,
    charges BLOB,
    key_value_pairs TEXT,
    data TEXT,
    natoms INTEGER)""",
    """CREATE TABLE species (
    Z INTEGER,
    n INTEGER,
    id INTEGER,
    FOREIGN KEY (id) REFERENCES systems(id))""",
    """CREATE TABLE keys (
    key TEXT,
    id INTEGER,
    FOREIGN KEY (id) REFERENCES systems(id))""",
    """CREATE TABLE text_key_values (
    key TEXT,
    value TEXT,
    id INTEGER,
    FOREIGN KEY (id) REFERENCES systems(id))""",
    """CREATE TABLE number_key_values (
    key TEXT,
    value REAL,
    id INTEGER,
    FOREIGN KEY (id) REFERENCES systems(id))""",
    """CREATE TABLE information (
    name TEXT,
    value TEXT)""",
    """INSERT INTO information VALUES ('version', '{0}')""".format(VERSION)]

index_statements = [
    'CREATE INDEX unique_id_index ON systems(unique_id)',
    'CREATE INDEX ctime_index ON systems(ctime)',
    'CREATE INDEX username_index ON systems(username)',
    'CREATE INDEX calculator_index ON systems(calculator)',
    'CREATE INDEX species_index ON species(Z)',
    'CREATE INDEX key_index ON keys(key)',
    'CREATE INDEX text_index ON text_key_values(key)',
    'CREATE INDEX number_index ON number_key_values(key)']

all_tables = ['systems', 'species', 'keys',
              'text_key_values', 'number_key_values']


class SQLite3Database(Database):
    initialized = False
    _allow_reading_old_format = False
    default = 'NULL'  # used for autoincrement id
    connection = None
    version = None
    
    def _connect(self):
        return sqlite3.connect(self.filename, timeout=600)

    def __enter__(self):
        self.connection = self._connect()
        return self
        
    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is None:
            self.connection.commit()
        else:
            self.connection.rollback()
        self.connection.close()
        self.connection = None
        
    def _initialize(self, con):
        if self.initialized:
            return

        cur = con.execute(
            'SELECT COUNT(*) FROM sqlite_master WHERE name="systems"')

        if cur.fetchone()[0] == 0:
            for statement in init_statements:
                con.execute(statement)
            if self.create_indices:
                for statement in index_statements:
                    con.execute(statement)
            con.commit()
            self.version = VERSION
        else:
            cur = con.execute(
                'SELECT COUNT(*) FROM sqlite_master WHERE name="user_index"')
            if cur.fetchone()[0] == 1:
                # Old version with "user" instead of "username" column
                self.version = 1
            else:
                try:
                    cur = con.execute(
                        'SELECT value FROM information WHERE name="version"')
                except sqlite3.OperationalError:
                    self.version = 2
                else:
                    self.version = int(cur.fetchone()[0])
                    
        self.initialized = True
                
    def _write(self, atoms, key_value_pairs, data):
        Database._write(self, atoms, key_value_pairs, data)
        
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
                
        id = None
        
        if isinstance(atoms, dict):
            dct = atoms
            unique_id = dct['unique_id']
            cur.execute('SELECT id FROM systems WHERE unique_id=?',
                        (unique_id,))
            rows = cur.fetchall()
            if rows:
                id = rows[0][0]
                self._delete(cur, [id], ['keys', 'text_key_values',
                                         'number_key_values'])
            dct['mtime'] = now()
        else:
            dct = self.collect_data(atoms)

        if 'constraints' in dct:
            constraints = encode(dct['constraints'])
        else:
            constraints = None
            
        numbers = dct.get('numbers')
        
        row = (dct['unique_id'],
               dct['ctime'],
               dct['mtime'],
               dct['user'],
               blob(numbers),
               blob(dct.get('positions')),
               blob(dct.get('cell')),
               int(np.dot(dct.get('pbc'), [1, 2, 4])),
               blob(dct.get('initial_magmoms')),
               blob(dct.get('initial_charges')),
               blob(dct.get('masses')),
               blob(dct.get('tags')),
               blob(dct.get('momenta')),
               constraints)

        if 'calculator' in dct:
            row += (dct['calculator'],
                    encode(dct['calculator_parameters']))
        else:
            row += (None, None)

        magmom = dct.get('magmom')
        if magmom is not None:
            # magmom can be one or three numbers (non-collinear case)
            magmom = np.array(magmom)
        row += (dct.get('energy'),
                dct.get('free_energy'),
                blob(dct.get('forces')),
                blob(dct.get('stress')),
                blob(dct.get('dipole')),
                blob(dct.get('magmoms')),
                blob(magmom),
                blob(dct.get('charges')),
                encode(key_value_pairs),
                encode(data),
                len(numbers))

        if id is None:
            q = self.default + ', ' + ', '.join('?' * len(row))
            cur.execute('INSERT INTO systems VALUES ({0})'.format(q),
                        row)
        else:
            q = ', '.join(line.split()[0].lstrip() + '=?'
                          for line in init_statements[0].splitlines()[2:])
            cur.execute('UPDATE systems SET {0} WHERE id=?'.format(q),
                        row + (id,))
        
        if id is None:
            id = self.get_last_id(cur)
            
            if len(numbers) > 0:
                count = np.bincount(numbers)
                unique_numbers = count.nonzero()[0]
                species = [(int(Z), int(count[Z]), id) for Z in unique_numbers]
                cur.executemany('INSERT INTO species VALUES (?, ?, ?)',
                                species)

        text_key_values = []
        number_key_values = []
        for key, value in key_value_pairs.items():
            if isinstance(value, (float, int)):
                number_key_values.append([key, float(value), id])
            else:
                assert isinstance(value, (str, unicode))
                text_key_values.append([key, value, id])
 
        cur.executemany('INSERT INTO text_key_values VALUES (?, ?, ?)',
                        text_key_values)
        cur.executemany('INSERT INTO number_key_values VALUES (?, ?, ?)',
                        number_key_values)
        cur.executemany('INSERT INTO keys VALUES (?, ?)',
                        [(key, id) for key in key_value_pairs])

        if self.connection is None:
            con.commit()
            con.close()
        return id
        
    def get_last_id(self, cur):
        cur.execute('SELECT seq FROM sqlite_sequence WHERE name="systems"')
        id = cur.fetchone()[0]
        return id
        
    def _get_dict(self, id):
        con = self._connect()
        c = con.cursor()
        if id is None:
            c.execute('SELECT COUNT(*) FROM systems')
            assert c.fetchone()[0] == 1
            c.execute('SELECT * FROM systems')
        else:
            c.execute('SELECT * FROM systems WHERE id=?', (id,))
        row = c.fetchone()
        return self.row_to_dict(row)

    def row_to_dict(self, row):
        if len(row) != 28:
            row = self._old2new(row)

        dct = {'id': row[0],
               'unique_id': row[1],
               'ctime': row[2],
               'mtime': row[3],
               'user': row[4],
               'numbers': deblob(row[5], np.int32),
               'positions': deblob(row[6], shape=(-1, 3)),
               'cell': deblob(row[7], shape=(3, 3)),
               'pbc': (row[8] & np.array([1, 2, 4])).astype(bool)}
        if row[9] is not None:
            dct['initial_magmoms'] = deblob(row[9])
        if row[10] is not None:
            dct['initial_charges'] = deblob(row[10])
        if row[11] is not None:
            dct['masses'] = deblob(row[11])
        if row[12] is not None:
            dct['tags'] = deblob(row[12], np.int32)
        if row[13] is not None:
            dct['momenta'] = deblob(row[13], shape=(-1, 3))
        if row[14] is not None:
            dct['constraints'] = decode(row[14])
        if row[15] is not None:
            dct['calculator'] = row[15]
            dct['calculator_parameters'] = decode(row[16])
        if row[17] is not None:
            dct['energy'] = row[17]
        if row[18] is not None:
            dct['free_energy'] = row[18]
        if row[19] is not None:
            dct['forces'] = deblob(row[19], shape=(-1, 3))
        if row[20] is not None:
            dct['stress'] = deblob(row[20])
        if row[21] is not None:
            dct['dipole'] = deblob(row[21])
        if row[22] is not None:
            dct['magmoms'] = deblob(row[22])
        if row[23] is not None:
            dct['magmom'] = deblob(row[23])[0]
        if row[24] is not None:
            dct['charges'] = deblob(row[24])

        for key, column in zip(['key_value_pairs', 'data'], row[25:27]):
            x = decode(column)
            if x:
                dct[key] = x
                
        return dct

    def _old2new(self, row):
        if not self._allow_reading_old_format:
            raise IOError('Please convert to new format. ' +
                          'Use: python -m ase.db.convert ' + self.filename)
        if len(row) == 26:
            extra = decode(row[25])
            return row[:-1] + (encode(extra['keywords']),
                               encode(extra['key_value_pairs']),
                               encode(extra['data']),
                               42)
        else:
            keywords = decode(row[-4])
            kvp = decode(row[-3])
            kvp.update(dict((keyword, 1) for keyword in keywords))
            return row[:-4] + (encode(kvp),) + row[-2:]
        
    def create_select_statement(self, keys, cmps, what='systems.*'):
        tables = ['systems']
        where = []
        args = []
        for n, key in enumerate(keys):
            tables.append('keys AS keys{0}'.format(n))
            where.append('systems.id=keys{0}.id AND keys{0}.key=?'.format(n))
            args.append(key)

        # Special handling of "H=0" and "H<2" type of selections:
        bad = {}
        for key, op, value in cmps:
            if isinstance(key, int):
                bad[key] = bad.get(key, True) and ops[op](0, value)
                
        nspecies = 0
        ntext = 0
        nnumber = 0
        for key, op, value in cmps:
            if key in ['id', 'energy', 'magmom', 'ctime', 'user',
                       'calculator', 'natoms']:
                if key == 'user' and self.version == 2:
                    key = 'username'
                where.append('systems.{0}{1}?'.format(key, op))
                args.append(value)
            elif isinstance(key, int):
                if bad[key]:
                    where.append(
                        'NOT EXISTS (SELECT id FROM species WHERE\n' +
                        '  species.id=systems.id AND species.Z==? AND ' +
                        'species.n{0}?)'.format(invop[op]))
                    args += [key, value]
                else:
                    tables.append('species AS specie{0}'.format(nspecies))
                    where.append(('systems.id=specie{0}.id AND ' +
                                  'specie{0}.Z=? AND ' +
                                  'specie{0}.n{1}?').format(nspecies, op))
                    args += [key, value]
                    nspecies += 1
            elif isinstance(value, (str, unicode)):
                tables.append('text_key_values AS text{0}'.format(ntext))
                where.append(('systems.id=text{0}.id AND ' +
                              'text{0}.key=? AND ' +
                              'text{0}.value{1}?').format(ntext, op))
                args += [key, value]
                ntext += 1
            else:
                tables.append('number_key_values AS number{0}'.format(nnumber))
                where.append(('systems.id=number{0}.id AND ' +
                              'number{0}.key=? AND ' +
                              'number{0}.value{1}?').format(nnumber, op))
                args += [key, float(value)]
                nnumber += 1
                
        sql = 'SELECT {0} FROM\n  '.format(what) + ', '.join(tables)
        if where:
            sql += '\n  WHERE\n  ' + ' AND\n  '.join(where)
        return sql, args
        
    def _select(self, keys, cmps, explain=False, verbosity=0,
                limit=None, offset=0):
        sql, args = self.create_select_statement(keys, cmps)
        
        if explain:
            sql = 'EXPLAIN QUERY PLAN ' + sql
            
        if limit:
            sql += '\nLIMIT {0}'.format(limit)

        if offset:
            sql += '\nOFFSET {0}'.format(offset)
            
        if verbosity == 2:
            print(sql, args)

        con = self._connect()
        self._initialize(con)
        cur = con.cursor()
        cur.execute(sql, args)
        if explain:
            for row in cur.fetchall():
                yield {'explain': row}
        else:
            for row in cur.fetchall():
                yield self.row_to_dict(row)
                    
    @parallel
    @lock
    def count(self, selection=None, **kwargs):
        keys, cmps = self.parse_selection(selection, **kwargs)
        sql, args = self.create_select_statement(keys, cmps, 'COUNT(*)')
        con = self._connect()
        self._initialize(con)
        cur = con.cursor()
        cur.execute(sql, args)
        return cur.fetchone()[0]
        
    def _update(self, ids, delete_keys, add_key_value_pairs):
        """Update row(s).
        
        ids: int or list of int
            ID's of rows to update.
        delete_keys: list of str
            Keys to remove.
        add_key_value_pairs: dict
            Key-value pairs to add.
            
        Returns number of key-value pairs added and keys removed.
        """
        
        dcts = [self._get_dict(id) for id in ids]
        m = 0
        n = 0
        with self:
            for dct in dcts:
                kvp = dct.get('key_value_pairs', {})
                n += len(kvp)
                for key in delete_keys:
                    kvp.pop(key, None)
                n -= len(kvp)
                m -= len(kvp)
                kvp.update(add_key_value_pairs)
                m += len(kvp)
                self._write(dct, kvp, data=dct.get('data', {}))
        return m, n

    @parallel
    @lock
    def delete(self, ids):
        con = self._connect()
        self._delete(con.cursor(), ids)
        con.commit()
        con.close()

    def _delete(self, cur, ids, tables=None):
        tables = tables or all_tables[::-1]
        for table in tables:
            cur.executemany('DELETE FROM {0} WHERE id=?'.format(table),
                            ((id,) for id in ids))


def blob(array):
    """Convert array to blob/buffer object."""

    if array is None:
        return None
    if array.dtype == np.int64:
        array = array.astype(np.int32)
    return buffer(np.ascontiguousarray(array))


def deblob(buf, dtype=float, shape=None):
    """Convert blob/buffer object to ndarray of correct dtype and shape.

    (without creating an extra view)."""

    if buf is None:
        return None
    if len(buf) == 0:
        array = np.zeros(0, dtype)
    else:
        array = np.frombuffer(buf, dtype)
    if shape is not None:
        array.shape = shape
    return array


if __name__ == '__main__':
    import sys
    from ase.db import connect
    con = connect(sys.argv[1])
    con._initialize(con._connect())
    print('Version:', con.version)
