from __future__ import print_function
import sys
import numpy as np

from ase.units import _hbar, _c, _e, Hartree
from gpaw.version import version
from gpaw.utilities.folder import Folder


def spectrum(exlist=None,
             filename=None,
             emin=None,
             emax=None,
             de=None,
             energyunit='eV',
             folding='Gauss',
             width=0.08,  # Gauss/Lorentz width
             comment=None,
             form='r'
             ):
    """Write out a folded spectrum.

    Parameters:
    =============== ===================================================
    ``exlist``      ExcitationList
    ``filename``    File name for the output file, STDOUT if not given
    ``emin``        min. energy, set to cover all energies if not given
    ``emax``        max. energy, set to cover all energies if not given
    ``de``          energy spacing
    ``energyunit``  Energy unit 'eV' or 'nm', default 'eV'
    ``folding``     Gauss (default) or Lorentz
    ``width``       folding width in terms of the chosen energyunit
    =============== ===================================================
    all energies in [eV]
    """

    # output
    out = sys.stdout
    if filename is not None:
        out = open(filename, 'w')
    if comment:
        print('#', comment, file=out)

    print('# Photoabsorption spectrum from linear response TD-DFT', file=out)
    print('# GPAW version:', version, file=out)
    if folding is not None:  # fold the spectrum
        print('# %s folded, width=%g [%s]' % (folding, width,
                                              energyunit), file=out)
    if form == 'r':
        out.write('# length form')
    else:
        assert(form == 'v')
        out.write('# velocity form')
    print('# om [%s]     osz          osz x       osz y       osz z'
          % energyunit, file=out)

    x = []
    y = []
    for ex in exlist:
        x.append(ex.get_energy() * Hartree)
        y.append(ex.get_oscillator_strength(form))

    if energyunit == 'nm':
        # transform to experimentally used wavelength [nm]
        x = 1.e+9 * 2 * np.pi * _hbar * _c / _e / np.array(x)
        y = np.array(y)
    elif energyunit != 'eV':
        raise RuntimeError('currently only eV and nm are supported')

    energies, values = Folder(width, folding).fold(x, y, de, emin, emax)
    for e, val in zip(energies, values):
        print('%10.5f %12.7e %12.7e %11.7e %11.7e' %
             (e, val[0], val[1], val[2], val[3]), file=out)

    if filename is not None:
        out.close()


def rotatory_spectrum(exlist=None,
                      filename=None,
                      emin=None,
                      emax=None,
                      de=None,
                      energyunit='eV',
                      folding='Gauss',
                      width=0.08,  # Gauss/Lorentz width
                      comment=None
                      ):
    """Write out a folded rotatory spectrum.

    See spectrum() for explanation of the parameters.
    """

    # output
    out = sys.stdout
    if filename is not None:
        out = open(filename, 'w')
    if comment:
        print('#', comment, file=out)

    print('# Rotatory spectrum from linear response TD-DFT', file=out)
    print('# GPAW version:', version, file=out)
    if folding is not None:  # fold the spectrum
        print('# %s folded, width=%g [%s]' % (folding, width,
                                              energyunit), file=out)
    print('# om [%s]     R [cgs]'
          % energyunit, file=out)

    x = []
    y = []
    for ex in exlist:
        x.append(ex.get_energy() * Hartree)
        y.append(ex.get_rotatory_strength())

    if energyunit == 'nm':
        # transform to experimentally used wavelength [nm]
        x = 1.e+9 * 2 * np.pi * _hbar * _c / _e / np.array(x)
        y = np.array(y)
    elif energyunit != 'eV':
        raise RuntimeError('currently only eV and nm are supported')

    energies, values = Folder(width, folding).fold(x, y, de, emin, emax)
    for e, val in zip(energies, values):
        print('%10.5f %12.7e' %
             (e, val), file=out)

    if filename is not None:
        out.close()


class Writer(Folder):

    def __init__(self, folding=None, width=0.08,  # Gauss/Lorentz width
                 ):
        self.folding = folding
        Folder.__init__(self, width, folding)

    def write(self, filename=None,
              emin=None, emax=None, de=None,
              comment=None):

        out = sys.stdout
        if filename is not None:
            out = open(filename, 'w')

        print('#', self.title, file=out)
        print('# GPAW version:', version, file=out)
        if comment:
            print('#', comment, file=out)
        if self.folding is not None:
            print('# %s folded, width=%g [eV]' % (self.folding,
                                                  self.width), file=out)
        print('#', self.fields, file=out)

        energies, values = self.fold(self.energies, self.values,
                                     de, emin, emax)
        for e, val in zip(energies, values):
            string = '%10.5f' % e
            for vf in val:
                string += ' %12.7e' % vf
            print(string, file=out)

        if filename is not None:
            out.close()
