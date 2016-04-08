import sys
import numpy as np
from numpy.linalg.linalg import LinAlgError
from ase.units import kJ
import ase.db
from ase.utils import prnt
from ase.test.tasks.dcdft import DeltaCodesDFTCollection as Collection
from ase.test.tasks.dcdft import FullEquationOfState as EquationOfState

collection = Collection()

db = sys.argv[1]

c = ase.db.connect(db)


def analyse(c, collection):

    A = []
    for name in collection.names:
        ve = []  # volume, energy pairs
        for d in c.select(name=name):
            try:
                ve.append((abs(np.linalg.det(d.cell)), d.energy))
            except AttributeError:
                ve.append((np.nan, np.nan))

        # sort according to volume
        ves = sorted(ve, key=lambda x: x[0])

        # EOS
        eos = EquationOfState([t[0] for t in ves],
                              [t[1] for t in ves])
        try:
            v, e, B0, B1, R = eos.fit()
        except (ValueError, TypeError, LinAlgError):
            (v, e, B0, B1, R) = (np.nan, np.nan, np.nan, np.nan, np.nan)
        e = e / len(collection[name])
        v = v / len(collection[name])
        B0 = B0 / kJ * 1.0e24  # GPa
        A.append((e, v, B0, B1, R))
    return np.array(A).T

E, V, B0, B1, R = analyse(c, collection)

fd = open(db + '_raw.txt', 'w')
for name, e, v, b0, b1, r, in zip(collection.names, E, V, B0, B1, R):
    if not np.isnan(e):
        prnt('%2s %8.4f %8.4f %8.4f' % (name, v, b0, b1), file=fd)
fd = open(db + '_raw.csv', 'w')
for name, e, v, b0, b1, r, in zip(collection.names, E, V, B0, B1, R):
    if not np.isnan(e):
        prnt('%s, %8.4f, %8.4f, %8.4f' % (name, v, b0, b1), file=fd)
