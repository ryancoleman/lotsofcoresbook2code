from __future__ import print_function
import sys
from optparse import OptionParser

def build_parser():
    description = 'Generate LCAO basis sets for the specified elements.'

    parser = OptionParser(usage='%prog [options] [elements]',
                          version='%prog 0.1', description=description)
    parser.add_option('-n', '--name', default=None, metavar='<name>',
                      help='name of generated basis files')
    parser.add_option('-t', '--type', default='dzp', metavar='<type>',
                      help='type of basis.  For example: sz, dzp, qztp, '+
                      '4z3p.  [default: %default]')
    parser.add_option('-E', '--energy-shift', metavar='<energy>', type='float',
                      default=.1,
                      help='use given energy shift to determine cutoff')
    parser.add_option('-T', '--tail-norm', metavar='<norm>', type='string',
                      default='0.16,0.3,0.6', dest='tailnorm',
                      help='use the given fractions to define the split'+
                      '-valence cutoffs.  Default: [%default]')
    parser.add_option('-f', '--xcfunctional', default='PBE', metavar='<XC>',
                      help='Exchange-Correlation functional '
                      '[default: %default]')
    parser.add_option('-g', '--non-relativistic-guess', action='store_true',
                      help='Run non-scalar relativistic AE calculation for '
                      'initial guess')
    parser.add_option('--rcut-max', type='float', default=16.,
                      metavar='<rcut>',
                      help='max cutoff for confined atomic orbitals.  This '+
                      'option has no effect on orbitals with smaller cutoff '+
                      '[default/Bohr: %default]')
    parser.add_option('--rcut-pol-rel', type='float', default=1.0,
                      metavar='<rcut>',
                      help='polarization function cutoff relative to largest'+
                      ' single-zeta cutoff [default: %default]')
    parser.add_option('--rchar-pol-rel', type='float', default=None,
                      metavar='<rchar>',
                      help='characteristic radius of Gaussian when not using '+
                      'interpolation scheme, relative to rcut')
    parser.add_option('--vconf-amplitude', type='float', default=12.,
                      metavar='<alpha>',
                      help='set proportionality constant of smooth '
                      'confinement potential [default: %default]')
    parser.add_option('--vconf-rstart-rel', type='float', default=.6,
                      metavar='<ri/rc>',
                      help='set inner cutoff for smooth confinement potential '
                      'relative to hard cutoff [default: %default]')
    parser.add_option('--vconf-sharp-confinement', action='store_true',
                      help='use sharp rather than smooth confinement '
                      'potential')
    parser.add_option('--debug', action='store_true',
                      help='use gpaw-DEBUG mode')
    parser.add_option('--lpol', type=int, default=None,
                      help='angular momentum quantum number '
                      'of polarization function.  '
                      'Default behaviour is to take the lowest l which is not '
                      'among the valence states.')
    parser.add_option('--jvalues',
                      help='explicitly specify which states to include.  '
                      'Numbering corresponds to generator\'s valence state '
                      'ordering.  '
                      'For example: 0,1,2.')

    return parser


bad_density_warning = """\
Bad initial electron density guess!  Try rerunning the basis generator
with the '-g' parameter to run a separate non-scalar relativistic
all-electron calculation and use its resulting density as an initial
guess."""

very_bad_density_warning = """\
Could not generate non-scalar relativistic electron density guess,
or non-scalar relativistic guess was not good enough for the scalar
relativistic calculation.  You probably have to use the Python interface
to the basis generator in gpaw.atom.basis directly and choose very
smart parameters."""


def main():
    parser = build_parser()
    opts, args = parser.parse_args()
    #if __name__ == '__main__' and len(args) == 0:
    #    parser.print_help()
    #    sys.exit(0)

    from gpaw.atom.basis import BasisMaker
    from gpaw import ConvergenceError
    from gpaw.basis_data import parse_basis_name

    zetacount, polcount, basistype = parse_basis_name(opts.type)
    referencefiles = [None] * len(args)
    reference_atom_indices = [None] * len(args)
    if polcount > 0:
        symbols = []
        for i, arg in enumerate(args):
            # Parse argument as <symbol>:<reference-file>:<nuclear index>.
            symbol_and_file = arg.split(':')
            symbol = symbol_and_file[0]
            symbols.append(symbol)
            if len(symbol_and_file) > 1:
                referencefiles[i] = symbol_and_file[1] # filename
            if len(symbol_and_file) == 3:
                reference_atom_indices[i] = int(symbol_and_file[2])
    else:
        symbols = args

    if opts.name is not None:
        name = '%s.%s' % (opts.name, basistype)
    else:
        name = basistype

    for symbol, referencefile, referenceindex in zip(symbols, referencefiles,
                                                     reference_atom_indices):
        try:
            bm = BasisMaker(symbol, name, gtxt=None,
                            non_relativistic_guess=opts.non_relativistic_guess,
                            xc=opts.xcfunctional)
        except ConvergenceError:
            if opts.non_relativistic_guess:
                print(very_bad_density_warning, file=sys.stderr)
                import traceback
                traceback.print_exc()
            else:
                print(bad_density_warning, file=sys.stderr)
            continue

        tailnorm = [float(norm) for norm in opts.tailnorm.split(',')]
        vconf_args = None
        if not opts.vconf_sharp_confinement:
            vconf_args = opts.vconf_amplitude, opts.vconf_rstart_rel

        jvalues = None
        if opts.jvalues:
            jvalues = [int(j) for j in opts.jvalues.split(',')]

        basis = bm.generate(zetacount, polcount,
                            tailnorm=tailnorm,
                            energysplit=opts.energy_shift,
                            referencefile=referencefile,
                            referenceindex=referenceindex,
                            rcutpol_rel=opts.rcut_pol_rel,
                            rcutmax=opts.rcut_max,
                            rcharpol_rel=opts.rchar_pol_rel,
                            vconf_args=vconf_args,
                            l_pol=opts.lpol,
                            jvalues=jvalues)
        basis.write_xml()
