from __future__ import print_function
from ase.io import write
from gpaw import restart

basename = 'CO'

# load nc binary file and get calculator
atoms, calc = restart(basename + '.gpw')

# write atomic positions to xyz-file
write(basename + '.xyz', atoms)

# loop over all wfs and write their cube files
nbands = calc.get_number_of_bands()

for band in range(nbands):
    wf = calc.get_pseudo_wave_function(band=band)
    fname=basename + '_' + '%d' % (band) + '.plt'
    print('writing wf', band, 'to file', fname)
    write(fname, atoms, data=wf)
