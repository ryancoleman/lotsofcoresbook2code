from ase import Atoms
from gpaw import GPAW
from ase.visualize import view


class ObsTest:
    def __init__(self):
        self.steps = []

    def run(self):
        self.steps.append(calc.iter)
        
    def converge(self):
        calc.scf.converged = True
            

H = Atoms('H', positions=[(1.5, 1.5, 1.5)])
H.set_cell((3, 3, 3))

calc = GPAW(convergence={'density': -1}, maxiter=11)
            
AllCall = ObsTest()
EvenCall = ObsTest()
OneCall = ObsTest()
FinalCall = ObsTest()

calc.attach(AllCall.run, 1)
calc.attach(EvenCall.run, 2)
calc.attach(OneCall.run, -7)
calc.attach(OneCall.converge, -7)
calc.attach(FinalCall.run, 0)

H.set_calculator(calc)
H.get_potential_energy()
    
assert AllCall.steps == [1, 2, 3, 4, 5, 6, 7]
assert EvenCall.steps == [2, 4, 6, 7]
assert OneCall.steps == [7]
assert FinalCall.steps == [7]
