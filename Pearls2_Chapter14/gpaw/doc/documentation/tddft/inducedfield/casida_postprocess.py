from gpaw import GPAW
from gpaw.lrtddft import LrTDDFT, photoabsorption_spectrum
from gpaw.inducedfield.inducedfield_lrtddft import LrTDDFTInducedField

# Load LrTDDFT object
lr = LrTDDFT('na2_lr.dat.gz')

# Calculate photoabsorption spectrum as usual
folding = 'Gauss'
width = 0.1
e_min = 0.0
e_max = 4.0
photoabsorption_spectrum(lr, 'na2_casida_spectrum.dat',
                         folding=folding, width=width,
                         e_min=e_min, e_max=e_max, delta_e=1e-2)

# Load GPAW object
calc = GPAW('na2_gs_casida.gpw')
print calc.wfs.rank_a

# Calculate induced field
frequencies = [1.0, 2.08]     # Frequencies of interest in eV
folding = 'Gauss'             # Folding function
width = 0.1                   # Line width for folding in eV
kickdir = 0                   # Kick field direction 0, 1, 2 for x, y, z
ind = LrTDDFTInducedField(paw=calc,
                          lr=lr,
                          frequencies=frequencies,
                          folding=folding,
                          width=width,
                          kickdir=kickdir)
ind.calculate_induced_field(gridrefinement=2, from_density='comp')
ind.write('na2_casida_field.ind', mode='field')
    
