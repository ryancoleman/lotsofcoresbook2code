#!/usr/bin/env python

import numpy as np

from ase.io import write
from gpaw.io.tar import Reader
from gpaw.tddft import TDDFT
from gpaw.tddft.fourier import DensityFourierTransform
from gpaw.tddft.units import aufrequency_to_eV, autime_to_attosec

# -------------------------------------------------------------------

if __name__ == '__main__':
    import os
    import sys

    def printOptions():
        scriptname = sys.argv[0].rsplit('/', 1)[-1]
        print('Usage:')
        print('    %s FTDFILE GPWFILE' % scriptname)
        print('')
        print('Arguments:')
        print('        FTDFILE    Fourier transformed density tar-file.')
        print('        GPWFILE    GPAW calculation tar-file (optional).')
        print('')

    try:
        assert len(sys.argv) == 3, 'Incorrect number of arguments.'

        ftd_filename = sys.argv[1]
        assert ftd_filename.endswith('.ftd'), 'Invalid FTD tarfile.'
        assert os.path.isfile(ftd_filename), 'FTD tarfile not found.'
        prefix = ftd_filename.rsplit('.ftd', 1)[0]

        tar = Reader(ftd_filename)
        try:
            timestep = tar['TimeStep']
            sigma = tar['Width']
        except KeyError:
            timestep = 1
            sigma = None
        omega_w = tar.get('Frequency')
        gamma_w = tar.get('PhaseFactor')
        Fnt_wsG = tar.get('FourierTransform')
        Ant_sG = tar.get('Average')
        atoms = None
        del tar

        gpw_filename = sys.argv[2]
        assert gpw_filename.endswith('.gpw'), 'Invalid GPW tarfile.'
        assert os.path.isfile(gpw_filename), 'GPW tarfile not found.'

        calc = TDDFT(gpw_filename, txt=None)
        obs = DensityFourierTransform(timestep * autime_to_attosec,
                                      omega_w * aufrequency_to_eV,
                                      (sigma is not None and sigma \
                                       * aufrequency_to_eV or None))
        obs.initialize(calc)
        atoms = calc.get_atoms()
        del calc
        obs.read(ftd_filename, idiotproof=False)
        try:
            sys.stdout.write('Select grid refinement [1*/2]: ')
            gdref = int(sys.stdin.readline().strip())
        except:
            gdref = 1
        getall = slice(None) #hack to obtain all frequencies/spins
        Fnt_wsG = obs.get_fourier_transform(getall, getall, gdref)
        Ant_sG = obs.get_average(getall, gdref)
        del obs

        # Save modulus and phase as .cube files for all frequencies/spins
        for w, Fnt_sG in enumerate(Fnt_wsG):
            for s, Fnt_G in enumerate(Fnt_sG):
                filename = '%s_Fnt_w%d_s%d_mod.cube' % (prefix,w,s)
                print('Saving %s (omega=%5.2f eV)...' \
                    % (filename, omega_w[w]*aufrequency_to_eV))
                write(filename, atoms, data=np.abs(Fnt_G))

                filename = '%s_Fnt_w%d_s%d_arg.cube' % (prefix,w,s)
                print('Saving %s (omega=%5.2f eV)...' \
                    % (filename, omega_w[w]*aufrequency_to_eV))
                write(filename, atoms, data=np.arctan2(Fnt_G.imag, Fnt_G.real))

        # Save mean density as .cube file for each spin
        for s, Ant_G in enumerate(Ant_sG):
            filename = '%s_Ant_s%d.cube' % (prefix,s)
            print('Saving %s...' % filename)
            write(filename, atoms, data=Ant_G)

    except AssertionError, e:
        printOptions()
        print('ERROR:', e)
        exit(-1)

