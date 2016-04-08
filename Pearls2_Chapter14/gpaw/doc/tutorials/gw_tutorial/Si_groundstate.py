from ase.lattice import bulk
from gpaw import GPAW, FermiDirac
from gpaw.wavefunctions.pw import PW

a = 5.431
atoms = bulk('Si', 'diamond', a=a)

calc = GPAW(
            mode=PW(200),                  # energy cutoff for plane wave basis (in eV)
            kpts=(3,3,3),
            xc='LDA',
            eigensolver='cg',
            occupations=FermiDirac(0.001),
            txt='Si_groundstate.txt'
           )

atoms.set_calculator(calc)
atoms.get_potential_energy()

calc.diagonalize_full_hamiltonian()       # determine all bands
calc.write('Si_groundstate.gpw','all')    # write out wavefunctions
