from ase.dft.wannier import Wannier
from gpaw import GPAW
from gpaw.mpi import world
calc = GPAW('Si-PBE.gpw', txt=None,
            parallel={'domain': world.size})
calc.wfs.ibz2bz(calc.atoms)
initial = [[(0.125, 0.125, 0.125), 0, 1.5],
           [(0.125, 0.625, 0.125), 0, 1.5],
           [(0.125, 0.125, 0.625), 0, 1.5],
           [(0.625, 0.125, 0.125), 0, 1.5]]
w = Wannier(4, calc,
            #fixedenergy=0.0,
            nbands=8,
            verbose=1,
            initialwannier=initial)
w.localize()
w.save('rotations.pckl')
