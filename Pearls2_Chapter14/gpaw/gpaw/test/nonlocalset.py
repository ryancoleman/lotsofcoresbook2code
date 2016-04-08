from ase import Atoms
from gpaw import GPAW
from gpaw.test import equal
from gpaw.poisson import PoissonSolver

atoms = Atoms('HF',positions=[(0.0,0.0,0.0),(1.0,0.0,0.0)])
atoms.set_pbc((True,True,True))
atoms.set_cell((2.0,2.0,2.0))

def MGGA_fail():
    calc = GPAW(xc='TPSS', eigensolver='cg',kpts=(1,2,1), convergence={'density':1e-8})
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    calc.set(kpts=(4,1,1))
    return atoms.get_potential_energy()

def MGGA_work():
    calc = GPAW(xc='TPSS', eigensolver='cg', kpts=(4,1,1), convergence={'density':1e-8})

    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

def GLLBSC_fail():
    calc = GPAW(xc='GLLBSC', eigensolver='cg', kpts=(1,2,1), convergence={'density':1e-8})

    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    calc.set(kpts=(4,1,1))
    return atoms.get_potential_energy()

def GLLBSC_work():
    calc = GPAW(xc='GLLBSC', eigensolver='cg', kpts=(4,1,1), convergence={'density':1e-8})

    atoms.set_calculator(calc)
    return atoms.get_potential_energy()

a = GLLBSC_fail()
b = GLLBSC_work()
equal(a,b, 1e-5)
c = MGGA_fail()
d = MGGA_work()
equal(c,d, 1e-5)

