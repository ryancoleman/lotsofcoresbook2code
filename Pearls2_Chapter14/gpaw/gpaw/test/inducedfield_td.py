from ase import Atoms
from gpaw import GPAW
from gpaw.tddft import TDDFT
from gpaw.inducedfield.inducedfield_tddft import TDDFTInducedField
import numpy as np

# Na2 cluster
atoms = Atoms(symbols='Na2', 
              positions=[(0, 0, 0), (3.0, 0, 0)],
              pbc=False)
atoms.center(vacuum=3.0)

# Standard ground state calculation
calc = GPAW(nbands=2, h=0.6, setups={'Na': '1'})
atoms.set_calculator(calc)
energy = atoms.get_potential_energy()
calc.write('na2_gs.gpw', mode='all')

# Standard time-propagation initialization
time_step = 10.0
iterations = 60
kick_strength = [1.0e-3, 1.0e-3, 0.0]
td_calc = TDDFT('na2_gs.gpw')
td_calc.absorption_kick(kick_strength=kick_strength)

# Create and attach InducedField object
frequencies = [1.0, 2.08]     # Frequencies of interest in eV
folding = 'Gauss'             # Folding function
width = 0.1                   # Line width for folding in eV
ind = TDDFTInducedField(paw=td_calc,
                        frequencies=frequencies,
                        folding=folding,
                        width=width,
                        restart_file='na2_td.ind')

# Propagate as usual
td_calc.propagate(time_step, iterations/2, 'na2_td_dm.dat', 'na2_td.gpw')

# Save TDDFT and InducedField objects
td_calc.write('na2_td.gpw', mode='all')
ind.write('na2_td.ind')

ind.paw = None

# Restart and continue
td_calc = TDDFT('na2_td.gpw')

# Load and attach InducedField object
ind = TDDFTInducedField(filename='na2_td.ind',
                        paw=td_calc,
                        restart_file='na2_td.ind')
    
# Continue propagation as usual
td_calc.propagate(time_step, iterations/2, 'na2_td_dm.dat', 'na2_td.gpw')

# Calculate induced electric field
ind.calculate_induced_field(gridrefinement=2, from_density='comp')

# Test
from gpaw.test import equal
tol  = 0.0001
val1 = ind.fieldgd.integrate(ind.Ffe_wg[0])
val2 = ind.fieldgd.integrate(np.abs(ind.Fef_wvg[0][0]))
val3 = ind.fieldgd.integrate(np.abs(ind.Fef_wvg[0][1]))
val4 = ind.fieldgd.integrate(np.abs(ind.Fef_wvg[0][2]))
val5 = ind.fieldgd.integrate(ind.Ffe_wg[1])
val6 = ind.fieldgd.integrate(np.abs(ind.Fef_wvg[1][0]))
val7 = ind.fieldgd.integrate(np.abs(ind.Fef_wvg[1][1]))
val8 = ind.fieldgd.integrate(np.abs(ind.Fef_wvg[1][2]))

equal(val1, 3307.77279745, tol)
equal(val2, 2.24614158834, tol)
equal(val3, 2.20381741663, tol)
equal(val4, 1.69244172329, tol)
equal(val5, 3305.78925228, tol)
equal(val6, 2.09432584636, tol)
equal(val7, 2.09849354306, tol)
equal(val8, 1.59727445257, tol)
    
ind.paw = None
