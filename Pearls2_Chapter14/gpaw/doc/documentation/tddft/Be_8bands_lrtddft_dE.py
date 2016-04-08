from gpaw import GPAW
from gpaw.lrtddft import LrTDDFT

c = GPAW('Be_gs_8bands.gpw')

dE = 10 # maximal Kohn-Sham transition energy to consider in eV
lr = LrTDDFT(c, xc='LDA', energy_range=dE)
lr.write('lr_dE.dat.gz')
