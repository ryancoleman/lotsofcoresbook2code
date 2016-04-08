from gpaw.tddft import TDDFT
from gpaw.inducedfield.inducedfield_tddft import TDDFTInducedField

# Load TDDFT object
td_calc = TDDFT('na2_td.gpw')

# Load and attach InducedField object
ind = TDDFTInducedField(filename='na2_td.ind',
                        paw=td_calc,
                        restart_file='na2_td.ind')

# Continue propagation as usual
time_step = 20.0
iterations = 250
td_calc.propagate(time_step, iterations, 'na2_td_dm.dat', 'na2_td.gpw')

# Save TDDFT and InducedField objects
td_calc.write('na2_td.gpw', mode='all')
ind.write('na2_td.ind')
