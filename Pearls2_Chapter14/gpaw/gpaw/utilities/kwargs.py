from __future__ import print_function
import os
import sys
    
from optparse import OptionParser, OptionGroup
from gpaw.parameters import InputParameters
from gpaw.poisson import PoissonSolver
from gpaw.mixer import Mixer, MixerSum, MixerDif
from gpaw.occupations import FermiDirac, ZeroKelvin
from gpaw.eigensolvers.rmm_diis import RMM_DIIS
from gpaw.eigensolvers.cg import CG
from gpaw.eigensolvers.davidson import Davidson
from gpaw.lcao.eigensolver import LCAO


def build_parser():
    usage = '%prog [OPTIONS] [FILE...]'
    description = ('Print representation of GPAW input parameters to stdout '
                   'or file.')
    parser = OptionParser(usage=usage, description=description)
    g = OptionGroup(parser, 'General')
    g.add_option('--complete', action='store_true', default=False,
                 help='print complete set of input parameters')
    g.add_option('--pop', metavar='KWARGS', default='',
                 help='comma-separated keyword arguments to pop')
    g.add_option('--output', metavar='FILE',
                 help='write to this file')
    g.add_option('--update', metavar='FILE',
                 help='update file in-place')
    parser.add_option_group(g)
    return parser


def kwargs2str(**kwargs):
    tokens = []
    start = 'dict('
    indent = ' ' * len(start)
    keys = kwargs.keys()
    keys.sort()
    for key in keys:
        tokens.append('%s=%s' % (key, repr(kwargs[key])))
    string = (',\n%s' % indent).join(tokens)
    return ''.join([start, string, ')'])


def str2kwargs(string):
    from numpy import array  # eval may need this
    kwargs = eval(string)
    assert isinstance(kwargs, dict)
    return kwargs


def append_to_optiongroup(parameters, opts):
    for key, value in parameters.items():
        opts.add_option('--%s' % key, default=repr(value), type=str,
                        help='default=%default')


def load(filename):
    string = open(filename).read()
    parameters = str2kwargs(string)
    return parameters


def populate_parser(parser, defaultparameters):
    opts = OptionGroup(parser, 'GPAW parameters')
    append_to_optiongroup(defaultparameters, opts)
    parser.add_option_group(opts)
    return parser


def main(argv):
    defaults = InputParameters()
    parser = build_parser()
    populate_parser(parser, defaults)

    opts, args = parser.parse_args(argv)

    if opts.update is not None and opts.output is not None:
        parser.error('Cannot both update and write to new file.')
    
    if opts.update:
        if not os.path.isfile(opts.update):
            parser.error('No such file %s' % opts.update)
        outfilename = opts.update
        args = [outfilename] + args
    elif opts.output:
        outfilename = opts.output
    else:
        outfilename = None
    
    parameters = {}
    if opts.complete:
        parameters.update(defaults)
    
    for arg in args:
        loaded_parameters = load(arg)
        parameters.update(loaded_parameters)
        # We have to use the newly loaded info (target file)
        # to get new defaults!
        #
        # Rather ugly but hopefully it works
        # We can probably avoid this somehow, think about it in the future
        # (can one have different call formats like e.g. the 'du' man-page,
        # and somehow detect which one is relevant?)
        #parser2 = build_parser()
        #populate_parser(parser2, loaded_parameters)
        #opts, args2 = parser2.parse_args(argv)
    
    cli_args = vars(opts)
    cli_parameters = {}
    keys = defaults.keys()
    for key in defaults:
        obj = eval(cli_args[key])
        if defaults[key] != obj:
            cli_parameters[key] = obj
    parameters.update(cli_parameters)

    # Remove keys which are not changed from defaults
    if not opts.complete:
        # must call keys() since dict cannot change size during iteration
        for key in parameters.keys(): 
            if defaults[key] == parameters[key]:
                del parameters[key]

    if opts.pop:
        for kwarg in opts.pop.split(','):
            parameters.pop(kwarg)
    
    if outfilename is None:
        out = sys.stdout
    else:
        out = open(outfilename, 'w')
    print(kwargs2str(**parameters), file=out)
