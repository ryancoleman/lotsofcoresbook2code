from ase.io import write
from gpaw import GPAW, Mixer
from ase.structure import molecule

CO = molecule('CO')
CO.set_cell((6., 6., 6.))
CO.center()

calc = GPAW(h=0.2,
            nbands=8,
            mixer=Mixer(beta=0.1, nmaxold=5, weight=50.0),
            txt='CO.txt')

CO.set_calculator(calc)
CO.get_potential_energy()

# Write wave functions to gpw file
calc.write('CO.gpw', mode='all')

# Generate cube-files of the orbitals.
for n in range(calc.get_number_of_bands()):
    wf = calc.get_pseudo_wave_function(band=n)
    write('CO%d.cube' % n, CO, data=wf)
