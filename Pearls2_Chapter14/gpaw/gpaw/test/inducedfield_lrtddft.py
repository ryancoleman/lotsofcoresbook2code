from ase import Atoms
from gpaw import GPAW
from gpaw.lrtddft import LrTDDFT, photoabsorption_spectrum
from gpaw.inducedfield.inducedfield_lrtddft import LrTDDFTInducedField
import numpy as np

# 1) Ground state calculation with empty states
atoms = Atoms(symbols='Na2', 
              positions=[(0, 0, 0), (3.0, 0, 0)],
              pbc=False)
atoms.center(vacuum=3.0)

calc = GPAW(nbands=20, h=0.6, setups={'Na': '1'})
atoms.set_calculator(calc)
energy = atoms.get_potential_energy()
calc.write('na2_gs_casida.gpw', mode='all')

# 2) Casida calculation
calc = GPAW('na2_gs_casida.gpw')
istart = 0
jend = 20
lr = LrTDDFT(calc, xc='LDA', istart=istart, jend=jend)
lr.diagonalize()
lr.write('na2_lr.dat.gz')

# 3) Calculate induced field
frequencies = [1.0, 2.08] # Frequencies of interest in eV
folding = 'Gauss'         # Folding function
width = 0.1               # Line width for folding in eV
kickdir = 0               # Kick field direction 0, 1, 2 for x, y, z
ind = LrTDDFTInducedField(paw=calc, lr=lr, frequencies=frequencies, folding=folding, width=width, kickdir=kickdir)
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
equal(val1, 3175.76177761, tol)
equal(val2, 1700.43442519, tol)
equal(val3, 1187.26249225, tol)
equal(val4, 1187.26249225, tol)
equal(val5, 10956.9813196, tol)
equal(val6, 6574.58868754, tol)
equal(val7, 4589.74440108, tol)
equal(val8, 4589.74440108, tol)

