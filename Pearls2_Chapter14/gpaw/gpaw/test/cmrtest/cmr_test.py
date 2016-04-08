import os
from ase import Atom, Atoms
from gpaw import GPAW
from gpaw.test import equal

import warnings
# cmr calls all available methods in ase.atoms detected by the module inspect.
# Therefore also deprecated methods are called - and we choose to silence those warnings.
warnings.filterwarnings('ignore', 'ase.atoms.*deprecated',)

import cmr

# from cmr.tools.log import Log
# cmr.logger.set_message_selection(Log.MSG_TYPE_ALL)


a = 4.05
d = a / 2 ** 0.5
bulk = Atoms([Atom('Al', (0, 0, 0)),
              Atom('Al', (0.5, 0.5, 0.5))],
             pbc=True)
bulk.set_cell((d, d, a), scale_atoms=True)
h = 0.18
calc = GPAW(h=h,
            nbands=2 * 8,
            kpts=(2, 2, 2),
            convergence={'energy': 1e-5})
bulk.set_calculator(calc)
e0 = bulk.get_potential_energy()

for ext in [".db", ".cmr"]:
    calc.write("test1"+ext)
    assert os.path.exists("test1"+ext)
    calc.write("test2"+ext, cmr_params={"value":1, "keywords":["a", "b"]})
    assert os.path.exists("test2"+ext)
    data = cmr.read("test2"+ext)
    assert data["value"] == 1
    assert len(data["db_keywords"]) == 2

# test the restart ability
if 0:
    calc = GPAW("test.db")
    calc.get_total_energy()
    calc = GPAW("test2.cmr")
    calc.get_total_energy()

