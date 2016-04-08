from __future__ import absolute_import, print_function
import os
import datetime
import json

import numpy as np

from ase.parallel import world
from ase.db.core import Database, ops, parallel, lock, now
from ase.db.core import reserved_keys


class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, datetime.datetime):
            return {'__datetime__': obj.isoformat()}
        if hasattr(obj, 'todict'):
            return obj.todict()
        return json.JSONEncoder.default(self, obj)


encode = MyEncoder().encode


def object_hook(dct):
    if '__datetime__' in dct:
        return datetime.datetime.strptime(dct['__datetime__'],
                                          '%Y-%m-%dT%H:%M:%S.%f')
    return dct


mydecode = json.JSONDecoder(object_hook=object_hook).decode


def intkey(key):
    try:
        return int(key)
    except ValueError:
        return key
    
    
def numpyfy(obj):
    if isinstance(obj, dict):
        return dict((intkey(key), numpyfy(value))
                    for key, value in obj.items())
    if isinstance(obj, list) and len(obj) > 0:
        try:
            a = np.array(obj)
        except ValueError:
            obj = [numpyfy(value) for value in obj]
        else:
            if a.dtype in [bool, int, float]:
                obj = a
    return obj


def decode(txt):
    return numpyfy(mydecode(txt))

    
def read_json(name):
    if isinstance(name, str):
        fd = open(name, 'r')
    else:
        fd = name
    dct = decode(fd.read())
    fd.close()
    return dct


class JSONDatabase(Database):
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_value, tb):
        pass
        
    def _write(self, atoms, key_value_pairs, data):
        Database._write(self, atoms, key_value_pairs, data)
        
        bigdct = {}
        ids = []
        nextid = 1

        if isinstance(self.filename, str) and os.path.isfile(self.filename):
            try:
                bigdct, ids, nextid = self._read_json()
            except (SyntaxError, ValueError):
                pass

        if isinstance(atoms, dict):
            dct = dict((key, atoms[key])
                       for key in reserved_keys
                       if key in atoms and key != 'id')
            unique_id = dct['unique_id']
            for id in ids:
                if bigdct[id]['unique_id'] == unique_id:
                    break
            else:
                id = None
            dct['mtime'] = now()
        else:
            dct = self.collect_data(atoms)
            id = None

        for key, value in [('key_value_pairs', key_value_pairs),
                           ('data', data)]:
            if value:
                dct[key] = value
            else:
                dct.pop(key, None)
        
        if id is None:
            id = nextid
            ids.append(id)
            nextid += 1
            
        bigdct[id] = dct
        self._write_json(bigdct, ids, nextid)
        return id
        
    def _read_json(self):
        bigdct = read_json(self.filename)
        ids = bigdct['ids']
        if not isinstance(ids, list):
            ids = ids.tolist()
        return bigdct, ids, bigdct['nextid']
        
    def _write_json(self, bigdct, ids, nextid):
        if world.rank > 0:
            return
            
        if isinstance(self.filename, str):
            fd = open(self.filename, 'w')
        else:
            fd = self.filename
        print('{', end='', file=fd)
        for id in ids:
            dct = bigdct[id]
            txt = ',\n '.join('"{0}": {1}'.format(key, encode(dct[key]))
                              for key in sorted(dct.keys()))
            print('"{0}": {{\n {1}}},'.format(id, txt), file=fd)
        print('"ids": {0},'.format(ids), file=fd)
        print('"nextid": {0}}}'.format(nextid), file=fd)

        if fd is not self.filename:
            fd.close()

    @parallel
    @lock
    def delete(self, ids):
        bigdct, myids, nextid = self._read_json()
        for id in ids:
            del bigdct[id]
            myids.remove(id)
        self._write_json(bigdct, myids, nextid)

    def _get_dict(self, id):
        bigdct, ids, nextid = self._read_json()
        if id is None:
            assert len(ids) == 1
            id = ids[0]
        dct = bigdct[id]
        dct['id'] = id
        return dct

    def _select(self, keys, cmps, explain=False, verbosity=0,
                limit=None, offset=0):
        if explain:
            yield {'explain': (0, 0, 0, 'scan table')}
            return
            
        try:
            bigdct, ids, nextid = self._read_json()
        except IOError:
            return
            
        if not limit:
            limit = -offset - 1
            
        cmps = [(key, ops[op], val) for key, op, val in cmps]
        n = 0
        for id in ids:
            if n - offset == limit:
                return
            dct = bigdct[id]
            for key in keys:
                if not ('key_value_pairs' in dct and
                        key in dct['key_value_pairs']):
                    break
            else:
                for key, op, val in cmps:
                    value = get_value(id, dct, key)
                    if not op(value, val):
                        break
                else:
                    if n >= offset:
                        dct['id'] = id
                        yield dct
                    n += 1

    def _update(self, ids, delete_keys, add_key_value_pairs):
        bigdct, myids, nextid = self._read_json()
        
        t = now()
        
        m = 0
        n = 0
        for id in ids:
            dct = bigdct[id]
            kvp = dct.get('key_value_pairs', {})
            n += len(kvp)
            for key in delete_keys:
                kvp.pop(key, None)
            n -= len(kvp)
            m -= len(kvp)
            kvp.update(add_key_value_pairs)
            m += len(kvp)
            if kvp:
                dct['key_value_pairs'] = kvp
            dct['mtime'] = t
            
        self._write_json(bigdct, myids, nextid)
        return m, n


def get_value(id, dct, key):
    pairs = dct.get('key_value_pairs')
    if pairs is None:
        value = None
    else:
        value = pairs.get(key)
    if value is not None:
        return value
    if key in ['energy', 'magmom', 'ctime', 'user', 'calculator']:
        return dct.get(key)
    if isinstance(key, int):
        return np.equal(dct['numbers'], key).sum()
    if key == 'natoms':
        return len(dct['numbers'])
    if key == 'id':
        return id
