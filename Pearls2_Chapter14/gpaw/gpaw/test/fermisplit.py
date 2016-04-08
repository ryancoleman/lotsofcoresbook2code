# this test should coverage the save and restore of
# fermi-levels when using fixmagmom:
#
# yes, fermi-level-splitting sounds a little bit strange
import numpy as np

from ase import Atoms
from gpaw import GPAW, FermiDirac, MixerSum
from gpaw.test import equal

modes = ['gpw']
try:
    import _gpaw_hdf5
    modes.append('hdf5')
except ImportError:
    pass

calc=GPAW(occupations=FermiDirac(width=0.1, fixmagmom=True),
          mixer=MixerSum(beta=0.05, nmaxold=3, weight=50.0),
          convergence={'energy':0.1, 'eigenstates': 1.5e-1,
                       'density': 1.5e-1})
atoms = Atoms('Cr', pbc=False)
atoms.center(vacuum=4)
mm = [1] * 1
mm[0] = 6.
atoms.set_initial_magnetic_moments(mm)
atoms.set_calculator(calc)
atoms.get_potential_energy()

ef1=calc.occupations.get_fermi_levels_mean()
efsplit1=calc.occupations.get_fermi_splitting()

ef3=calc.occupations.get_fermi_levels()
for mode in modes:
    calc.write('test.%s' % mode)

for mode in modes:
    # check number one: is the splitting value saved?
    readtest=GPAW("test.%s" % mode)
    ef2=readtest.occupations.get_fermi_levels_mean()
    efsplit2=readtest.occupations.get_fermi_splitting()
    
    # numpy arrays
    ef4=readtest.occupations.get_fermi_levels()

    # These values should be identic
    equal(ef1, ef2, 1e-9)
    equal(efsplit1,  efsplit2, 1e-9)
    equal(ef3.mean(), ef1, 1e-9)
    equal(ef3.mean(), ef2, 1e-9)
    equal(ef3.mean(), ef4.mean(), 1e-9)
    equal(ef3[0] - ef3[1], ef4[0] - ef4[1], 1e-9)
    equal(efsplit1, ef4[0] - ef4[1], 1e-9)
