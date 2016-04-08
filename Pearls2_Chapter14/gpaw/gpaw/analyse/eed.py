from __future__ import print_function
import sys

import numpy as np
from ase.units import Bohr, Hartree
from ase.parallel import paropen
from ase.data import vdw_radii

import _gpaw
from gpaw.io.fmf import FMF


class ExteriorElectronDensity:

    """Exterior electron density to describe MIES spectra.

    Simple approach to describe MIES spectra after
    Y. Harada et al., Chem. Rev. 97 (1997) 1897
    """

    def __init__(self, gd, atoms):
        """Find the grid points outside of the van der Waals radii
        of the atoms"""

        assert gd.orthogonal
        self.gd = gd

        n = len(atoms)
        atom_c = atoms.positions / Bohr
        vdWradius = np.empty((n))
        for a, atom in enumerate(atoms):
            vdWradius[a] = self.get_vdWradius(atom.number)

        # define the exterior region mask
        mask = gd.empty(dtype=int)
        _gpaw.eed_region(mask, atom_c, gd.beg_c, gd.end_c,
                         gd.h_cv.diagonal().copy(), vdWradius)
        self.mask = mask

    def get_weight(self, psit_G):
        """Get the weight of a wave function in the exterior region
        (outside of the van der Waals radius). The augmentation sphere
        is assumed to be smaller than the van der Waals radius and hence
        does not contribute."""

        # smooth part
        weigth = self.gd.integrate(np.where(self.mask == 1,
                                            (psit_G * psit_G.conj()).real,
                                            0.0))

        return weigth

    def get_vdWradius(self, Z):
        """Return van der Waals radius in Bohr"""
        r = vdw_radii[Z] / Bohr
        if np.isnan(r):
            msg = 'van der Waals radius for Z=' + str(Z) + ' not known!'
            raise RuntimeError(msg)
        else:
            return r

    def write_mies_weights(self, wfs, file=None):
        if file is None:
            file = 'eed_mies.dat'

        if isinstance(file, str):
            out = paropen(file, 'aw')
        else:
            out = file

        fmf = FMF(['exterior electron density weights after',
                   'Y. Harada et al., Chem. Rev. 97 (1997) 1897'])
        print(fmf.header(), end=' ', file=out)
        print(fmf.data(['band index: n',
                        'k-point index: k',
                        'spin index: s',
                        'k-point weight: weight',
                        'energy: energy [eV]',
                        'occupation number: occ',
                        'relative EED weight: eed_weight']), end=' ', file=out)

        print(
            '#; n   k s   weight      energy         occ  eed_weight', file=out)
        for kpt in wfs.kpt_u:
            for n in range(wfs.bd.nbands):
                print('%4d %3d %1d %8.5f  %10.5f  %10.5f  %10.5f' %
                     (n, kpt.k, kpt.s, kpt.weight,
                      kpt.eps_n[n] * Hartree,
                      kpt.f_n[n],
                      self.get_weight(kpt.psit_nG[n])
                      ), file=out)
                if hasattr(out, 'flush'):
                    out.flush()
