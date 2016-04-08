import os
import warnings
# cmr calls all available methods in ase.atoms detected by the module inspect.
# Therefore also deprecated methods are called - and we choose to silence those warnings.
warnings.filterwarnings('ignore', 'ase.atoms.*deprecated',)

from ase.test import NotAvailable

# if CMR_SETTINGS_FILE is missing, cmr raises simply
# Exception("CMR is not configured properly. Please create the settings file with cmr --create-settings.")
try:
    import cmr
except (Exception, ImportError):
    raise NotAvailable('CMR is required')

from ase.calculators.emt import EMT
from ase.io import read, write
from ase.structure import molecule

cmr_params = {"db_keywords":["O", "ase"], # keyword
              "molecule":"O2"} #field

m1 = molecule('O2')
m1.set_calculator(EMT())
e1 = m1.get_potential_energy()
write("O2.cmr", m1, cmr_params = cmr_params)

reread = read("O2.cmr")
e2 = reread.get_potential_energy()
assert abs(e1 - e2) < 1.e-6, str(e1) + ' ' + str(e2)

db_read = cmr.read("O2.cmr")
assert "O" in db_read["db_keywords"]
assert "ase" in db_read["db_keywords"]
assert db_read["molecule"] == "O2"

# clean
filename = "O2.cmr"
if os.path.exists(filename): os.unlink(filename)
