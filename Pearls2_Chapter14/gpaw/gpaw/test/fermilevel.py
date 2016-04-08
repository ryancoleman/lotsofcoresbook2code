import numpy as np

from ase import Atoms
from gpaw import GPAW, FermiDirac
from gpaw.test import equal

modes = ['gpw']
try:
    import _gpaw_hdf5
    modes.append('hdf5')
except ImportError:
    pass

calc = GPAW(nbands=1, occupations=FermiDirac(0.0))#, txt=None)
atoms = Atoms('He', pbc=True, calculator=calc)
atoms.center(vacuum=3)

e0 = atoms.get_potential_energy()
niter0 = calc.get_number_of_iterations()
try:
    calc.get_fermi_level()
except ValueError:
    pass # It *should* raise an error
else:
    raise RuntimeError, 'get_fermi_level should not be possible for width=0'
calc.set(nbands=3, convergence={'bands':2})
atoms.get_potential_energy()
homo, lumo = calc.get_homo_lumo()
equal(homo, -15.4473, 0.01)
equal(lumo,  -0.2566, 0.01)
for mode in modes:
    calc.write('test.%s' % mode)
    assert np.all(GPAW('test.%s' % mode, txt=None).get_homo_lumo() == (homo, lumo))
    ef = calc.get_fermi_level()
    equal(ef, -7.85196, 0.01)

calc.set(occupations=FermiDirac(0.1))
e1 = atoms.get_potential_energy()
niter1 = calc.get_number_of_iterations()
ef = calc.get_fermi_level()
equal(ef, -7.85196, 0.01)
for mode in modes:
    calc.write('test.%s' % mode)
    equal(GPAW('test.%s' % mode, txt=None).get_fermi_level(), ef, 1e-8)
