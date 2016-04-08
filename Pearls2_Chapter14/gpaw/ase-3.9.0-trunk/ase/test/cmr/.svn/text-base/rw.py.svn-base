import os
import warnings
# cmr calls all available methods in ase.atoms detected by the module inspect.
# Therefore also deprecated methods are called - and we choose to silence those warnings.
warnings.filterwarnings('ignore', 'ase.atoms.*deprecated',)

import numpy as np

def array_almost_equal(a1, a2, tol=np.finfo(type(1.0)).eps):
    """Replacement for old numpy.testing.utils.array_almost_equal."""
    return (np.abs(a1 - a2) < tol).all()

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

m1 = molecule('O2')
m1.center(2.0)

write("O2.cmr", images=m1)

m1.set_calculator(EMT())
e1 = m1.get_potential_energy()
f1 = m1.get_forces()

m2 = read("O2.cmr")

m2.set_calculator(EMT())
e2 = m2.get_potential_energy()
f2 = m1.get_forces()

# assume atoms definitions are the same if energy/forces are the same: can we do better?
assert abs(e1-e2) < 1.e-6, str(e1) + ' ' + str(e2)
assert array_almost_equal(f1, f2, tol=1.e-6)

# clean
filename = "O2.cmr"
if os.path.exists(filename): os.unlink(filename)

