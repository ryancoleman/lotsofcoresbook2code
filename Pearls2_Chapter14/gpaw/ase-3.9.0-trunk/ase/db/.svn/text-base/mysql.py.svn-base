import MySQLdb

from ase.db.sqlite import init_statements, index_statements
from ase.db.sqlite import all_tables, SQLite3Database


class Connection:
    def __init__(self, con):
        self.con = con

    def cursor(self):
        return Cursor(self.con.cursor())
    
    def commit(self):
        self.con.commit()

    def close(self):
        self.con.close()


class Cursor:
    def __init__(self, cur):
        self.cur = cur

    def fetchone(self):
        return self.cur.fetchone()

    def fetchall(self):
        return self.cur.fetchall()

    def execute(self, statement, *args):
        self.cur.execute(statement.replace('?', '%s'), *args)

    def executemany(self, statement, *args):
        self.cur.executemany(statement.replace('?', '%s'), *args)

    
class MySQLDatabase(SQLite3Database):
    default = 'DEFAULT'
    
    def _connect(self):
        con = MySQLdb.connect(db='mysql', user='ase', passwd='ase')
        return Connection(con)

    def _initialize(self, con):
        pass
    
    def get_last_id(self, cur):
        cur.execute('select max(id) from systems')
        id = cur.fetchone()[0]
        return int(id)


def reset():
    db = 'mysql'
    con = MySQLdb.connect(db=db, user='root')
    cur = con.cursor()

    cur.execute("use " + db)
    cur.execute("show tables like 'systems'")
    if cur.fetchone() is not None:
        # order matters for drop (drop 'systems' last)
        for t in ['information'] + all_tables[::-1]:  # MDTMP is information special?
            cur.execute("drop table %s.%s" % (db, t))
            # avoid "Commands out of sync; you can't run this command now"
            # by closing the cursor after each execute
            cur.close()
            cur = con.cursor()
        cur.execute("drop user 'ase'@'localhost'")
    cur.execute("create user 'ase'@'localhost' identified by 'ase'")
    con.commit()
    # MySQL can grant privileges only for one table at a time?
    for t in all_tables[:] + ['information']:   # MDTMP is information special?
        cur.execute("grant all privileges on %s.%s to 'ase'@'localhost'" % (db, t))
        cur.close()
        cur = con.cursor()
    con.commit()

    for sql in init_statements:
        for a, b in [('BLOB', 'BLOB'),
                     ('TEXT', 'VARCHAR(767)'),
                     ('data VARCHAR(767)', 'data TEXT'),
                     ('REAL', 'DOUBLE'),
                     ('KEY AUTOINCREMENT',
                      'KEY AUTO_INCREMENT')]:
            sql = sql.replace(a, b)
        cur.execute(sql)
        cur.close()
        cur = con.cursor()
    for sql in index_statements:
        cur.execute(sql)
        cur.close()
        cur = con.cursor()
    con.commit()


if __name__ == '__main__':
    # Debian
    # PYTHONPATH=/path/to/ase python -m ase.db.mysql
    # RHEL:
    # su -c "yum -y remove mariadb-server"
    # su -c "rm -rf /var/lib/mysql"
    # su -c "yum -y install mariadb-server"
    # su -c "systemctl start mariadb.service"
    # PYTHONPATH=/path/to/ase python -m ase.db.mysql
    reset()
