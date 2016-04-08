from __future__ import print_function
import os

from ase.units import Hartree

from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft.spectrum import spectrum, rotatory_spectrum


def check_convergence(lr,             # LrTDDFT object
                      dirname='conv',  # directory name to store the files
                      title=None,     # title for gnuplot
                      dn=None,        # steps to vary istart/jend
                      dE=None,        # steps to vary energy range
                      emin=None,
                      emax=None,
                      folding='Gauss',
                      istart=None,
                      jend=None,
                      width=0.03):     # Gauss/Lorentz width
    """Study the convergence of a LrTDDFT calcualtion by varying istart/jend.
    A gnuplot file will be created with the name 'dirname'/conv.gpl."""

    if istart is None:
        istart0 = lr.kss.istart
    else:
        if istart < lr.kss.istart:
            raise RuntimeError
        istart0 = istart
    if jend is None:
        jend0 = lr.kss.jend
    else:
        if jend > lr.kss.jend:
            raise RuntimeError
        jend0 = jend

    # create subdirectory for the files
    if not os.path.isdir(dirname):
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        else:
            raise RuntimeError('Can\'t create directory ' + dirname)

    def fname(filename):
        return dirname + '/' + filename

    fgpl = open(fname('conv.gpl'), 'w')
    print('set xlabel "omega [eV]"', file=fgpl)
    print('set ylabel "Folded osc. strength [1/eV]"', file=fgpl)
    if not emin:
        emin_gpl = '*'
    else:
        emin_gpl = str(emin)
    if not emax:
        emax_gpl = '*'
    else:
        emax_gpl = str(emax)
    print('set xrange [' + emin_gpl + ':' + emax_gpl + ']', file=fgpl)
    if title:
        print('set title "' + str(title) + '"', file=fgpl)

    # kss
    spectrum(lr.kss, fname('kss.dat'), width=width)
    spectrum(lr.kss, fname('ksssticks.dat'), folding=None)

    # full
    lr.diagonalize(istart=istart0, jend=jend0)
    spectrum(lr, fname('full.dat'), width=width)
    spectrum(lr, fname('fullsticks.dat'), folding=None)
    print('plot "' + fname('full.dat') + '" t "full" w l lt 1, \\', file=fgpl)
    print('     "' + fname('fullsticks.dat') +
          '" u 1:($2*20) t "" w impulses lt 1, \\', file=fgpl)
    print('     "' + fname('kss.dat') + '" t "Kohn-Sham" w l lt 2', file=fgpl)
    print('pause -10', file=fgpl)

    if dn is None:
        dn = -istart0 + jend0
        dn = int(dn / 10.)

    if dE is None:
        # vary istart
        print('plot "' + fname('full.dat') + '" t "istart=' +
              str(istart0) + '" w l lt 1, \\', file=fgpl)
        for i in range(1, 4):
            istart = istart0 + i * dn
            lr.diagonalize(istart=istart, jend=jend0)
            fn = fname('istart' + str(istart) + '.dat')
            spectrum(lr, fn, width=width)
            print('    "' + fn + '" t "istart=' + str(istart)
                  + '" w l lt', i + 1, end=' ', file=fgpl)
            if i < 3:
                print(', \\', end=' ', file=fgpl)
            print(file=fgpl)
        print('pause -10', file=fgpl)

        # vary jend
        print('plot "' + fname('full.dat') + '" t "jend=' +
              str(jend0) + '" w l lt 1, \\', file=fgpl)
        for i in range(1, 4):
            jend = jend0 - i * dn
            lr.diagonalize(jend=jend, istart=istart0)
            fn = fname('jend' + str(jend) + '.dat')
            spectrum(lr, fn, width=width)
            print('    "' + fn + '" t "jend=' + str(jend)
                  + '" w l lt', i + 1, end=' ', file=fgpl)
            if i < 3:
                print(', \\', end=' ', file=fgpl)
            print(file=fgpl)
        print('pause -10', file=fgpl)
    else:
        # vary the enrgy range
        max_kss_energy = 0.
        for kss in lr.Om.kss:
            max_kss_energy = max(max_kss_energy, kss.get_energy() * Hartree)
        print('plot "' + fname('full.dat') + '" t "full"' +
              ' w l lt 1, \\', file=fgpl)
        for i in range(1, 4):
            max_energy = max_kss_energy - i * dE
            lr.diagonalize(energy_range=max_energy)
            fn = fname('max_energy' + str(max_energy) + '.dat')
            spectrum(lr, fn, width=width)
            print('    "' + fn + '" t "dE=-' + str(
                i * dE) + '" w l lt', end=' ', file=fgpl)
            print(i + 1, end=' ', file=fgpl)
            if i < 3:
                print(', \\', end=' ', file=fgpl)
            print(file=fgpl)
        print('pause -10', file=fgpl)

    # plot different directions
    print('plot "' + fname('full.dat') +
          '" u 1:3 t "x" w l lt 1, \\', file=fgpl)
    print('     "' + fname('full.dat') +
          '" u 1:4 t "y" w l lt 2, \\', file=fgpl)
    print('     "' + fname('full.dat') + '" u 1:5 t "z" w l lt 3', file=fgpl)
    print('pause -10', file=fgpl)

    # plot rotary strength
    if lr[0].magn is not None:
        print('set ylabel "Folded rot. strength [cgs/eV]"', file=fgpl)
        rotatory_spectrum(lr, fname('rotatory.dat'), width=width)
        rotatory_spectrum(lr, fname('rotatory_sticks.dat'), folding=None)
        print('plot "' + fname('rotatory.dat') +
              '" t "rotatory" w l lt 1, \\', file=fgpl)
        print('     "' + fname('rotatory_sticks.dat') +
              '" u 1:($2*20) t "" w impulses lt 1', file=fgpl)
        print('pause -10', file=fgpl)

    fgpl.close()
