# Emacs: treat this as -*- python -*-

import os
from optparse import OptionParser

def build_parser():
    usage = '%prog [OPTION] [BASIS]...'
    parser = OptionParser(usage=usage, version='%prog 1.0')
    parser.add_option('-f', '--files', action='store_true',
                      dest='actual_filenames',
                      help='Read from specified filenames rather than '
                      'searching GPAW setup directories')
    parser.add_option('-s', '--save-figs', action='store_true', dest='save',
                      help='Save figures to disk rather than showing plots')
    parser.add_option('-l', '--literal', action='store_true',
                      help='Do not pre-multiply wave functions by r in plots')
    parser.add_option('-n', '--normalize', action='store_true',
                      help='Plot normalized wave functions')
    parser.add_option('-x', '--ext', default='png',
                      help='Image format [default: %default]')
    return parser

def main():
    parser = build_parser()
    opts, files = parser.parse_args()

    import pylab
    from gpaw.basis_data import Basis, BasisPlotter

    plotter = BasisPlotter(premultiply=not opts.literal,
                           normalize=opts.normalize,
                           show=False,
                           save=opts.save,
                           ext=opts.ext)

    for path in files:
        dir, filename = os.path.split(path)
        splitfilename = filename.split('.')
        symbol = splitfilename[0]
        name = '.'.join(splitfilename[1:-1])
        if opts.actual_filenames:
            basis = Basis(symbol, name, False)
            basis.read_xml(path)
        else: # Search GPAW setup dirs
            basis = Basis(symbol, name)        
        plotter.plot(basis)

    if not opts.save:
        pylab.show()
