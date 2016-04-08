from ase import Atoms
from gpaw import GPAW, ConvergenceError
from gpaw.mixer import Mixer, MixerSum

# C2 singlet from http://toc.uni-muenster.de/GMTKN/GMTKN30/W4-08.html
m = Atoms(symbols='C2',
          positions=[
    [ 0.        ,  0.        , -0.62000006],
    [ 0.        ,  0.        ,  0.62000006]],
          )
m.center(vacuum=4.0)

calc = GPAW(h=0.18,
            xc='PBE',
            basis='dzp',
            maxiter=550,
            width=0.0,
            txt='C2_default.txt',
            )

m.set_calculator(calc)
try:
    m.get_potential_energy()
except ConvergenceError:
    pass

assert not calc.scf.converged

del calc

# converges (fortuitously)
# C2 really needs broken spin-symmetry to converge
calc = GPAW(h=0.18,
            xc='PBE',
            basis='dzp',
            maxiter=550,
            width=0.0,
            )
calc.set(
    mixer=MixerSum(0.02, 3),
    eigensolver='cg',
    spinpol=True,
    txt='C2_conv1.txt',
    )

m.set_calculator(calc)
try:
    e1 = m.get_potential_energy()
except ConvergenceError:
    e1 = None
    pass

del calc

# or with broken symmetry magnetic moments, to a different solution
m.set_initial_magnetic_moments([-0.5,0.5])
calc = GPAW(h=0.18,
            xc='PBE',
            basis='dzp',
            maxiter=550,
            width=0.0,
            )
calc.set(
    txt='C2_conv2.txt',
    )

m.set_calculator(calc)
e2 = m.get_potential_energy()

if e1 is not None:
    # Note that spin-symmetry broken solution gives a different energy!
    # Standard DFT is unable to treat such systems, similarly to the
    # famous H2 dissociation: dx.doi.org/10.1103/PhysRevLett.87.133004
    assert e2 - e1 > 0.15
