from ase.io import read
from gpaw import GPAW, Mixer, ConvergenceError
from gpaw.eigensolvers.rmm_diis import RMM_DIIS
from gpaw.mpi import size
from gpaw import use_mic

atoms = read('POSCAR')

if use_mic:
    txt = 'out_C60_Pb100_mic_p%d.txt' % size
else:
    txt = 'out_C60_Pb100_p%d.txt' % size

calc = GPAW(h=0.18, nbands=-180, width=0.2,
            kpts=(2,2,1), xc='PBE',
            maxiter=15,
            txt=txt, eigensolver=RMM_DIIS(niter=2),
            mixer=Mixer(0.1, 5, 100),
            # parallel={'sl_default': (4,4,64)}
           )
atoms.set_calculator(calc)
try:
    atoms.get_potential_energy()
except ConvergenceError:
    pass
