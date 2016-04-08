from ase import Atoms
from gpaw import GPAW, FermiDirac, Mixer
#from gpaw.xc.noncolinear import NonColinearLDA, NonColinearLCAOEigensolver, \
#     NonColinearMixer

h = Atoms('H', magmoms=[1])
h.center(vacuum=2)
xc = 'LDA'
c = GPAW(txt='c.txt',
         mode='lcao',
         basis='dz(dzp)',
         #setups='ncpp',
         h=0.25,
         xc=xc,
         #occupations=FermiDirac(0.01),
         mixer=Mixer(),
         #noncolinear=[(2,0,0)],
         )#eigensolver=NonColinearLCAOEigensolver())
c.set(nbands=1)
h.calc = c
h.get_potential_energy()
