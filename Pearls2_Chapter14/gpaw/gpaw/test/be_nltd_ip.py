import numpy as np
from ase import Atoms
from gpaw import GPAW
from gpaw.tddft import TDDFT
from gpaw.tddft.abc import LinearAbsorbingBoundary
from gpaw.tddft.laser import CWField

atoms = Atoms('Be',[(0,0,0)], pbc=False)
atoms.center(vacuum=6)
calc = GPAW(h=0.35)
atoms.set_calculator(calc)
atoms.get_potential_energy()

calc.write('be_gs.gpw', 'all')


td_calc = TDDFT('be_gs.gpw',
                td_potential = CWField(1e-3, 2.0*np.pi/50.0, 150.0))
td_calc.set_absorbing_boundary(LinearAbsorbingBoundary(5.0, 0.01,
                                                       atoms.positions.copy()))
td_calc.propagate(8.0, 5,
                  'be_nl_dmz_ipabs_1e-3.dat',
                  'be_nl_td.gpw')

td_rest = TDDFT('be_nl_td.gpw',
                td_potential = CWField(1e-3, 2.0*np.pi/50.0, 150.0))
td_rest.set_absorbing_boundary(LinearAbsorbingBoundary(5.0, 0.01,
                                                       atoms.positions.copy()))
td_rest.propagate(8.0, 5,
                  'be_nl_dmz_ipabs_1e-3.dat',
                  'be_nl_td.gpw')
