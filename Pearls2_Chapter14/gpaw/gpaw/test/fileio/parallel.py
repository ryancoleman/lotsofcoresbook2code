from gpaw import GPAW, ConvergenceError, restart
from gpaw.eigensolvers.rmm_diis_old import RMM_DIIS
from gpaw.mixer import MixerSum
from gpaw.test import equal
from ase.lattice import bulk

modes = ['gpw']
try:
    import _gpaw_hdf5
    modes.append('hdf5')
except ImportError:
    pass

# bulk Fe with k-point, band, and domain parallelization
a = 2.87
atoms = bulk('Fe', 'bcc', a=a)
atoms.set_initial_magnetic_moments([2.2,])
calc = GPAW(h=0.20,
            eigensolver=RMM_DIIS(),
            mixer=MixerSum(0.1,3),
            nbands=6,
            kpts=(4,4,4),
            parallel={'band' : 2, 'domain' : (2,1,1)},
            maxiter=4)
atoms.set_calculator(calc)
try:
    atoms.get_potential_energy()
except ConvergenceError:
    pass
for mode in modes:
    calc.write('tmp.%s' % mode, mode='all')

# Continue calculation for few iterations
for mode in modes:
    atoms, calc = restart('tmp.%s' % mode,
                          eigensolver=RMM_DIIS(),
                          mixer=MixerSum(0.1,3),
                          parallel={'band' : 2, 'domain' : (1,1,2)},
                          maxiter=4)
    try:
        atoms.get_potential_energy()
    except ConvergenceError:
        pass
    e = calc.hamiltonian.Etot
    equal(e, -0.372602008394, 0.000001)
