import ase.db
c = ase.db.connect('results.db')
failed = ', '.join(d.name for d in c.select(ok=0))
assert not failed, 'Failed to converge: ' + failed

