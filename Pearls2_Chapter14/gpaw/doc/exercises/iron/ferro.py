from ase import Atoms
from gpaw import GPAW, PW

a = 2.87
m = 2.2

fe = Atoms('Fe2',
           scaled_positions=[(0, 0, 0),
                             (0.5, 0.5, 0.5)],
           magmoms=[m, m],
           cell=(a, a, a),
           pbc=True)

calc = GPAW(mode=PW(350),
            kpts=(6, 6, 6),
            txt='ferro.txt')

fe.set_calculator(calc)
e = fe.get_potential_energy()
calc.write('ferro.gpw')
