from gpaw import GPAW
from gpaw.lrtddft import LrTDDFT

c = GPAW('Be_gs_8bands.gpw')

istart=0 # band index of the first occ. band to consider
jend=8  # band index of the last unocc. band to consider
lr = LrTDDFT(c, xc='LDA', istart=istart, jend=jend,
             nspins=2) # force the calculation of triplet excitations also
lr.write('lr.dat.gz')
