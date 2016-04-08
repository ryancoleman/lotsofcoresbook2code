from ase import Atoms

from gpaw import GPAW, FermiDirac

k = 7
kT = 0
h = 0.18

a = 5.4
b = a/2
atoms = Atoms('Si2',
              positions=([0, 0, 0],
                         [b/2, b/2, b/2]),
              cell=([0, b, b],
                    [b, 0, b],
                    [b, b, 0]),
              pbc=True)
   
calc = GPAW(kpts=(k, k, k),
            setups='ah',
            symmetry='off',
            occupations=FermiDirac(kT),
            h=h)

atoms.set_calculator(calc)
atoms.get_potential_energy()
calc.write('Si_AH.gpw', mode='all')
