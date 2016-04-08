import os
import sys

from ase.db import connect
from ase.db.sqlite import index_statements


def convert(name):
    con1 = connect(name, use_lock_file=False)
    con1._allow_reading_old_format = True
    newname = name[:-2] + 'new.db'
    with connect(newname, create_indices=False, use_lock_file=False) as con2:
        for dct in con1.select():
            kvp = dct.get('key_value_pairs', {})
            con2.write(dct, data=dct.get('data'), **kvp)
        
    c = con2._connect()
    for statement in index_statements:
        c.execute(statement)
    c.commit()

    os.rename(name, name[:-2] + 'old.db')
    os.rename(newname, name)
    
        
if __name__ == '__main__':
    for name in sys.argv[1:]:
        convert(name)
