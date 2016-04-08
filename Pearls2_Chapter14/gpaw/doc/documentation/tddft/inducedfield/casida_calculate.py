from ase import Atoms
from gpaw import GPAW
from gpaw.lrtddft import LrTDDFT

# Na2 cluster
atoms = Atoms(symbols='Na2', 
              positions=[(0, 0, 0), (3.0, 0, 0)],
              pbc=False)
atoms.center(vacuum=6.0)

# Standard ground state calculation with empty states
calc = GPAW(nbands=100, h=0.4, setups={'Na': '1'})
atoms.set_calculator(calc)
energy = atoms.get_potential_energy()

calc.set(convergence={'bands' : 90}, 
         fixdensity=True,
         eigensolver='cg')
calc.calculate()

calc.write('na2_gs_casida.gpw', mode='all')

# Standard Casida calculation
calc = GPAW('na2_gs_casida.gpw')
istart = 0
jend = 90
lr = LrTDDFT(calc, xc='LDA', istart=istart, jend=jend)
lr.diagonalize()
lr.write('na2_lr.dat.gz')
