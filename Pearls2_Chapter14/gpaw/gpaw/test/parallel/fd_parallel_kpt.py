from __future__ import print_function
import sys

from ase import Atoms
from ase.utils import devnull

from gpaw import GPAW, FermiDirac
from gpaw import KohnShamConvergenceError
from gpaw.utilities import compiled_with_sl

from ase.structure import molecule

# Calculates energy and forces for various parallelizations

tolerance = 4e-5

parallel = dict()

basekwargs = dict(mode='fd',
                  eigensolver='rmm-diis',
                  maxiter=3,
                  #basis='dzp',
                  #nbands=18,
                  nbands=6,
                  kpts=(4,4,4), # 8 kpts in the IBZ
                  parallel=parallel)

Eref = None
Fref_av = None


def run(formula='H2O', vacuum=1.5, cell=None, pbc=1, **morekwargs):
    print(formula, parallel)
    system = molecule(formula)
    kwargs = dict(basekwargs)
    kwargs.update(morekwargs)
    calc = GPAW(**kwargs)
    system.set_calculator(calc)
    system.center(vacuum)
    if cell is None:
        system.center(vacuum)
    else:
        system.set_cell(cell)
    system.set_pbc(pbc)

    try:
        system.get_potential_energy()
    except KohnShamConvergenceError:
        pass

    E = calc.hamiltonian.Etot
    F_av = calc.forces.calculate(calc.wfs, calc.density,
                                 calc.hamiltonian)

    global Eref, Fref_av
    if Eref is None:
        Eref = E
        Fref_av = F_av

    eerr = abs(E - Eref)
    ferr = abs(F_av - Fref_av).max()

    if calc.wfs.world.rank == 0:
        print('Energy', E)
        print()
        print('Forces')
        print(F_av)
        print()
        print('Errs', eerr, ferr)

    if eerr > tolerance or ferr > tolerance:
        if calc.wfs.world.rank == 0:
            stderr = sys.stderr
        else:
            stderr = devnull
        if eerr > tolerance:
            print('Failed!', file=stderr)
            print('E = %f, Eref = %f' % (E, Eref), file=stderr)
            msg = 'Energy err larger than tolerance: %f' % eerr
        if ferr > tolerance:
            print('Failed!', file=stderr)
            print('Forces:', file=stderr)
            print(F_av, file=stderr)
            print(file=stderr)
            print('Ref forces:', file=stderr)
            print(Fref_av, file=stderr)
            print(file=stderr)
            msg = 'Force err larger than tolerance: %f' % ferr
        print(file=stderr)
        print('Args:', file=stderr)
        print(formula, vacuum, cell, pbc, morekwargs, file=stderr)
        print(parallel, file=stderr)
        raise AssertionError(msg)
        

# reference:
# kpt-parallelization = 8,
# state-parallelization = 1,
# domain-decomposition = (1,1,1)
run()

# kpt-parallelization = 2,
# state-parallelization = 2,
# domain-decomposition = (1,2,1)
parallel['band'] = 2
parallel['domain'] = (1, 2, 1)
run()

if compiled_with_sl():
    # kpt-parallelization = 2,
    # state-parallelization = 2,
    # domain-decomposition = (1,2,1)
    # with blacs
    parallel['sl_default'] = (2, 2, 2)
    run()

# perform spin polarization test
parallel = dict()

basekwargs = dict(mode='fd',
                  eigensolver='rmm-diis',
                  maxiter=3,
                  nbands=6,
                  kpts=(4,4,4), # 8 kpts in the IBZ
                  parallel=parallel)

Eref = None
Fref_av = None

OH_kwargs = dict(formula='NH2', vacuum=1.5, pbc=1, spinpol=1,
                 occupations=FermiDirac(width=0.1))

# reference:
# kpt-parallelization = 4,
# spin-polarization = 2,
# state-parallelization = 1,
# domain-decomposition = (1,1,1)
run(**OH_kwargs)

# kpt-parallelization = 2,
# spin-polarization = 2,
# domain-decomposition = (1, 2, 1)
parallel['domain'] = (1, 2, 1)
run(**OH_kwargs)

# kpt-parallelization = 2,
# spin-polarization = 2,
# state-parallelization = 2,
# domain-decomposition = (1, 1, 1)
del parallel['domain']
parallel['band'] = 2
run(**OH_kwargs) 

# do last test plus buffer_size keyword
parallel['buffer_size'] = 150
run(**OH_kwargs)

if compiled_with_sl():
    # kpt-parallelization = 2,
    # spin-polarization = 2,
    # state-parallelization = 2,
    # domain-decomposition = (1, 2, 1)
    # with blacs
    del parallel['buffer_size']
    parallel['domain'] = (1, 2, 1)
    parallel['sl_default'] = (2, 1, 2)
    run(**OH_kwargs)

    # kpt-parallelization = 2,
    # state-parallelization = 2,
    # domain-decomposition = (1, 2, 1)
    # with blacs
    parallel['sl_default'] = (2, 2, 2)
    run(**OH_kwargs)

    # do last test plus buffer_size keyword
    parallel['buffer_size'] = 150
    run(**OH_kwargs)

