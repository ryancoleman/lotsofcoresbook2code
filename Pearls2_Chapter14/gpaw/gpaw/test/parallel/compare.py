from ase import Atoms, Atom
from gpaw import GPAW
from gpaw.mpi import rank, size
a = 3.0
H = Atoms([Atom('H')],
          cell=(a, a, a),
          pbc=True,
          calculator=GPAW())
if size > 1:
    H.positions[0, 0] += 0.01 * rank
    try:
        e0 = H.get_potential_energy()
    except RuntimeError:
        pass
    else:
        assert False
