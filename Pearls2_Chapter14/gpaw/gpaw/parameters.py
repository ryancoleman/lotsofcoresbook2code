import numpy as np

from ase.units import Hartree, Bohr
from ase.dft.kpoints import monkhorst_pack

import gpaw
import gpaw.mpi as mpi
from gpaw.wavefunctions.pw import PW
from gpaw.occupations import FermiDirac
from gpaw.poisson import PoissonSolver, FFTPoissonSolver


def usesymm2symmetry(usesymm):
    if usesymm is None:
        return {'time_reversal': False, 'point_group': False}
    if usesymm:
        return {'time_reversal': True, 'point_group': True}
    return {'time_reversal': True, 'point_group': False}


class InputParameters(dict):
    def __init__(self, **kwargs):
        dict.__init__(self, [
            ('h', None),  # Angstrom
            ('xc', 'LDA'),
            ('gpts', None),
            ('kpts', [(0, 0, 0)]),
            ('lmax', 2),
            ('charge', 0),
            ('fixmom', False),  # don't use this
            ('nbands', None),
            ('setups', 'paw'),
            ('basis', {}),
            ('width', None),  # eV, don't use this
            ('occupations', None),
            ('spinpol', None),
            ('usesymm', 'default'),  # don't use this
            ('stencils', (3, 3)),
            ('fixdensity', False),
            ('mixer', None),
            ('txt', '-'),
            ('hund', False),
            ('random', False),
            ('dtype', None),
            ('filter', None),
            ('maxiter', 333),  # google it's spiritual meaning!
            ('parallel', {'kpt': None,
                          'domain': gpaw.parsize_domain,
                          'band': gpaw.parsize_bands,
                          'order': 'kdb',
                          'stridebands': False,
                          'sl_auto': False,
                          'sl_default': gpaw.sl_default,
                          'sl_diagonalize': gpaw.sl_diagonalize,
                          'sl_inverse_cholesky': gpaw.sl_inverse_cholesky,
                          'sl_lcao': gpaw.sl_lcao,
                          'sl_lrtddft': gpaw.sl_lrtddft,
                          'buffer_size': gpaw.buffer_size}),
            ('parsize', None),  # don't use this
            ('parsize_bands', None),  # don't use this
            ('parstride_bands', False),  # don't use this
            ('external', None),  # eV
            ('verbose', 0),
            ('eigensolver', None),
            ('poissonsolver', None),
            ('communicator', mpi.world),
            ('idiotproof', True),
            ('mode', 'fd'),
            ('convergence', {'energy': 0.0005,  # eV / electron
                             'density': 1.0e-4,
                             'eigenstates': 4.0e-8,  # eV^2
                             'bands': 'occupied',
                             'forces': None}),  # eV / Ang Max
            ('realspace', None),
            ('symmetry', {'point_group': True,
                          'time_reversal': True,
                          'symmorphic': True,
                          'tolerance': 1e-7})])
        dict.update(self, kwargs)

    def __repr__(self):
        dictrepr = dict.__repr__(self)
        repr = 'InputParameters(**%s)' % dictrepr
        return repr
    
    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        assert key in self
        self[key] = value

    def update(self, parameters):
        for key, value in parameters.items():
            assert key in self
        dict.update(self, parameters)

    def read(self, reader):
        """Read state from file."""

        r = reader

        version = r['version']
        
        assert version >= 0.3
    
        self.xc = r['XCFunctional']
        self.nbands = r.dimension('nbands')
        self.spinpol = (r.dimension('nspins') == 2)

        bzk_kc = r.get('BZKPoints', broadcast=True)
        if r.has_array('NBZKPoints'):
            self.kpts = r.get('NBZKPoints', broadcast=True)
            if r.has_array('MonkhorstPackOffset'):
                offset_c = r.get('MonkhorstPackOffset', broadcast=True)
                if offset_c.any():
                    self.kpts = monkhorst_pack(self.kpts) + offset_c
        else:
            self.kpts = bzk_kc

        if version < 4:
            self.symmetry = usesymm2symmetry(r['UseSymmetry'])
        else:
            self.symmetry = {'point_group': r['SymmetryOnSwitch'],
                             'symmorphic': r['SymmetrySymmorphicSwitch'],
                             'time_reversal': r['SymmetryTimeReversalSwitch'],
                             'tolerance': r['SymmetryToleranceCriterion']}

        try:
            self.basis = r['BasisSet']
        except KeyError:
            pass

        if version >= 2:
            try:
                h = r['GridSpacing']
            except KeyError:  # CMR can't handle None!
                h = None
            if h is not None:
                self.h = Bohr * h
            if r.has_array('GridPoints'):
                self.gpts = r.get('GridPoints')
        else:
            if version >= 0.9:
                h = r['GridSpacing']
            else:
                h = None

            gpts = ((r.dimension('ngptsx') + 1) // 2 * 2,
                    (r.dimension('ngptsy') + 1) // 2 * 2,
                    (r.dimension('ngptsz') + 1) // 2 * 2)

            if h is None:
                self.gpts = gpts
            else:
                self.h = Bohr * h

        self.lmax = r['MaximumAngularMomentum']
        self.setups = r['SetupTypes']
        self.fixdensity = r['FixDensity']
        if version <= 0.4:
            # Old version: XXX
            print(('# Warning: Reading old version 0.3/0.4 restart files ' +
                  'will be disabled some day in the future!'))
            self.convergence['eigenstates'] = r['Tolerance']
        else:
            nbtc = r['NumberOfBandsToConverge']
            if not isinstance(nbtc, (int, str)):
                # The string 'all' was eval'ed to the all() function!
                nbtc = 'all'
            if version < 5:
                force_crit = None
            else:
                force_crit = r['ForcesConvergenceCriterion']
                if force_crit is not None:
                    force_crit *= (Hartree / Bohr)
            self.convergence = {'density': r['DensityConvergenceCriterion'],
                                'energy':
                                r['EnergyConvergenceCriterion'] * Hartree,
                                'eigenstates':
                                r['EigenstatesConvergenceCriterion'],
                                'bands': nbtc,
                                'forces': force_crit}

            if version < 1:
                # Volume per grid-point:
                dv = (abs(np.linalg.det(r.get('UnitCell'))) /
                      (gpts[0] * gpts[1] * gpts[2]))
                self.convergence['eigenstates'] *= Hartree**2 * dv

            if version <= 0.6:
                mixer = 'Mixer'
                weight = r['MixMetric']
            elif version <= 0.7:
                mixer = r['MixClass']
                weight = r['MixWeight']
                metric = r['MixMetric']
                if metric is None:
                    weight = 1.0
            else:
                mixer = r['MixClass']
                weight = r['MixWeight']

            if mixer == 'Mixer':
                from gpaw.mixer import Mixer
            elif mixer == 'MixerSum':
                from gpaw.mixer import MixerSum as Mixer
            elif mixer == 'MixerSum2':
                from gpaw.mixer import MixerSum2 as Mixer
            elif mixer == 'MixerDif':
                from gpaw.mixer import MixerDif as Mixer
            elif mixer == 'DummyMixer':
                from gpaw.mixer import DummyMixer as Mixer
            else:
                Mixer = None

            if Mixer is None:
                self.mixer = None
            else:
                self.mixer = Mixer(r['MixBeta'], r['MixOld'], weight)
            
        if version == 0.3:
            # Old version: XXX
            print(('# Warning: Reading old version 0.3 restart files is ' +
                  'dangerous and will be disabled some day in the future!'))
            self.stencils = (2, 3)
            self.charge = 0.0
            fixmom = False
        else:
            self.stencils = (r['KohnShamStencil'],
                             r['InterpolationStencil'])
            if r['PoissonStencil'] == 999:
                self.poissonsolver = FFTPoissonSolver()
            else:
                self.poissonsolver = PoissonSolver(nn=r['PoissonStencil'])
            self.charge = r['Charge']
            fixmom = r['FixMagneticMoment']

        self.occupations = FermiDirac(r['FermiWidth'] * Hartree,
                                      fixmagmom=fixmom)

        try:
            self.mode = r['Mode']
        except KeyError:
            self.mode = 'fd'

        if self.mode == 'pw':
            self.mode = PW(ecut=r['PlaneWaveCutoff'] * Hartree)
            
        if len(bzk_kc) == 1 and not bzk_kc[0].any():
            # Gamma point only:
            if r['DataType'] == 'Complex':
                self.dtype = complex
