"""Module defining  ``Eigensolver`` classes."""

from gpaw.eigensolvers.rmm_diis import RMM_DIIS
from gpaw.eigensolvers.cg import CG
from gpaw.eigensolvers.davidson import Davidson
from gpaw.lcao.eigensolver import DirectLCAO
from gpaw import use_mic


def get_eigensolver(name, mode, convergence=None):
    """Create eigensolver object."""
    if name is None:
        if mode.name == 'lcao':
            name = 'lcao'
        else:
            name = 'dav'
    if isinstance(name, str):
        eigensolver = {'rmm-diis': RMM_DIIS,
                       'cg': CG,
                       'dav': Davidson,
                       'lcao': DirectLCAO
                       }[name]()
    else:
        eigensolver = name
    
    if use_mic and not isinstance(eigensolver, RMM_DIIS):
        raise NotImplementedError('MIC version works only with RMM_DIIS')
    if isinstance(eigensolver, CG):
        eigensolver.tolerance = convergence['eigenstates']

    assert isinstance(eigensolver, DirectLCAO) == (mode.name == 'lcao')

    return eigensolver
