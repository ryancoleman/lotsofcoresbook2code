from ase import Atom, io, optimize
from gpaw import GPAW, FermiDirac
from gpaw.cluster import Cluster
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft.excited_state import ExcitedState

box = 5.     # box dimension
h = 0.25     # grid spacing
width = 0.01 # Fermi width
nbands = 6   # bands in GS calculation
nconv = 4    # bands in GS calculation to converge
R = 2.99     # starting distance
iex = 1      # excited state index
d = 0.01     # step for numerical force evaluation
exc = 'LDA'  # xc for the linear response TDDFT kernel

s = Cluster([Atom('Na'), Atom('Na', [0, 0, R])])
s.minimal_box(box, h=h)

c = GPAW(h=h, nbands=nbands, eigensolver='cg',
         occupations=FermiDirac(width=width),
         setups={'Na': '1'},
         convergence={'bands':nconv})
c.calculate(s)
lr = LrTDDFT(c, xc=exc, eps=0.1, jend=nconv-1)

ex = ExcitedState(lr, iex, d=d)
s.set_calculator(ex)

ftraj='relax_ex' + str(iex)
ftraj += '_box' + str(box) + '_h' + str(h)
ftraj += '_d' + str(d) + '.traj'
traj = io.PickleTrajectory(ftraj, 'w', s)
dyn = optimize.FIRE(s)
dyn.attach(traj.write)
dyn.run(fmax=0.05)
