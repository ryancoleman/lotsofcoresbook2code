from __future__ import print_function
from optparse import OptionParser


def build_parser():
    parser = OptionParser(usage='%prog [options] [elements]',
                          version='%prog 0.1')
    parser.add_option('-f', '--xcfunctional', type='string', default='LDA',
                      help='Exchange-Correlation functional ' +
                      '(default value LDA)', metavar='<XC>')
    parser.add_option('-n', '--non-scalar-relativistic', action='store_true',
                      default=False,
                      help='Do *not* do a scalar-relativistic calculation.')
    parser.add_option('-x', '--exact-exchange', action='store_true',
                      default=False,
                      help='Calculate exact exchange integrals.')
    parser.add_option('-r', '--radius', type='string', default=None,
                      help='Cutoff radius or radii (comma separated).',
                      metavar='<rcut>')
    parser.add_option('-v', '--zero-potential', metavar='type,radius',
                      help='Type of zero-potential - type must be either ' +
                      '"poly" or "f".')
    parser.add_option('--filter', metavar='h,x',
                      help='Parameters used for Fourier-filtering and '
                      'projector functions and zero-potential. "h" is '
                      'the cutoff grid-spacing (in Bohr) and "x" is the ratio '
                      'between outer and inner radii.')
    parser.add_option('-l', '--logarithmic-derivatives', action='store_true',
                      help='Calculate logarithmic derivatives.')
    parser.add_option('-a', '--all-electron-only', action='store_true',
                      help='Skip generation of PAW setup.')
    parser.add_option('-e', '--extra-projectors', type='string', default=None,
                      help='Extra projectors. Use ";" to separate s, p and ' +
                      'd channels.  Examples: "0.0,1.0" for two extra ' +
                      's-type. "0.0;1.0" for extra s and p. ";1.0" for ' +
                      'extra p.', metavar='0.0;0.0,1.0;0.0')
    parser.add_option('-c', '--core', type='string', default=None,
                      help='Frozen core.  Examples: "[Ne]", "[Ar]3d".',
                      metavar='<core>')
    parser.add_option('--normconserving', type='string',
                      help='Examples: s, sp.')
    parser.add_option('--core-hole', metavar='state,occ',
                      help='Add core hole. Examples: "1s,0.5", "2p,1".')
    parser.add_option('--configuration', metavar='config',
                      help='Specify non-groundstate configuration. '
                      'Na+ ion: "Ne,3s0", O2- ion: "1s2,2s2,2p6" or ' +
                      '"He,2s2,2p6".')
    parser.add_option('--compensation-charge-radius', metavar='rcut',
                      type=float,
                      help='Cutoff radius for compensation charges.')
    parser.add_option('--name', type='string', metavar='<id>',
                      help='Name to use for setup file: <symbol>.<id>.<xc>.  '
                      'Default name is <symbol>.<xc>.')
    parser.add_option('--use-restart-file', action='store_true',
                      default=False,
                      help='Use restart file (should be avoided: '
                      'introduces dependency on the restart).')
    parser.add_option('-w', '--write-files', action='store_true',
                      help='Write wave functions and other things to files.')
    parser.add_option('-p', '--plot', action='store_true',
                      help='Show plot and generate reStructuredText.')
    parser.add_option('-g', '--points-per-node', metavar='<gpernode>',
                      type=int, default=150,
                      help='Number of radial grid points per node.')
    parser.add_option('--empty-states', type='string', default=None,
                      help='Add empty state(s).  Example: 5p.',
                      metavar='<states>')
    parser.add_option('--tf-coefficient', type='float', default=1,
                      help='Sets value of coefficient in Thomas-Fermi ' +
                      'calculations. Default is 1',
                      metavar='<tf-coefficient>')
    parser.add_option('--orbital-free', action='store_true',
                      help='Generates orbital-free Thomas-Fermi setup')
    return parser


def main():
    parser = build_parser()
    opt, args = parser.parse_args()

    import sys
    
    from gpaw.atom.generator import Generator
    from gpaw.atom.configurations import parameters, tf_parameters
    from gpaw.atom.all_electron import AllElectron
    from gpaw import ConvergenceError

    if args:
        atoms = args
    else:
        atoms = parameters.keys()

    bad_density_warning = """\
    Problem with initial electron density guess!  Try to run the program
    with the '-nw' option (non-scalar-relativistic calculation + write
    density) and then try again without the '-n' option (this will
    generate a good initial guess for the density)."""

    for symbol in atoms:
        scalarrel = not opt.non_scalar_relativistic

        corehole = None
        if opt.core_hole is not None:
            state, occ = opt.core_hole.split(',')
            # Translate corestate string ('1s') to n and l:
            ncorehole = int(state[0])
            lcorehole = 'spdf'.find(state[1])
            fcorehole = float(occ)
            corehole = (ncorehole, lcorehole, fcorehole)

        if opt.all_electron_only:
            a = AllElectron(symbol, opt.xcfunctional, scalarrel, corehole,
                            opt.configuration, not opt.write_files, '-',
                            opt.points_per_node,
                            opt.orbital_free, opt.tf_coefficient)
            try:
                a.run()
            except ConvergenceError:
                print(bad_density_warning, file=sys.stderr)
            continue
        g = Generator(symbol, opt.xcfunctional, scalarrel, corehole,
                      opt.configuration, not opt.write_files, '-',
                      opt.points_per_node, orbital_free=opt.orbital_free,
                      tf_coeff=opt.tf_coefficient)

        if opt.orbital_free:
            p = tf_parameters.get(symbol, {'rcut': 0.9})
        else:
            p = parameters.get(symbol, {})

        if opt.core is not None:
            p['core'] = opt.core

        if opt.radius is not None:
            p['rcut'] = [float(x) for x in opt.radius.split(',')]

        if opt.extra_projectors is not None:
            extra = {}
            if opt.extra_projectors != '':
                for l, x in enumerate(opt.extra_projectors.split(';')):
                    if x != '':
                        extra[l] = [float(y) for y in x.split(',')]
            p['extra'] = extra

        if opt.normconserving is not None:
            p['normconserving'] = opt.normconserving

        if opt.filter is not None:
            p['filter'] = [float(x) for x in opt.filter.split(',')]

        if opt.compensation_charge_radius is not None:
            p['rcutcomp'] = opt.compensation_charge_radius

        if opt.zero_potential is not None:
            vbar = opt.zero_potential.split(',')
            p['vbar'] = (vbar[0], float(vbar[1]))

        if opt.empty_states is not None:
            p['empty_states'] = opt.empty_states

        try:
            g.run(logderiv=opt.logarithmic_derivatives,
                  exx=opt.exact_exchange, name=opt.name,
                  use_restart_file=opt.use_restart_file,
                  **p)
        except ConvergenceError:
            print(bad_density_warning, file=sys.stderr)
        except RuntimeError, m:
            if len(m.__str__()) == 0:
                raise
            print(m)

        if opt.plot:
            from gpaw.atom.analyse_setup import analyse
            analyse(g, show=True)
