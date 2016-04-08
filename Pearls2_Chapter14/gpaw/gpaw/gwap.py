"""GPAW Without Any Python (GWAP)."""
from __future__ import print_function
import os
import sys
import inspect
import optparse
import textwrap


functions = {'xc': 'gpaw.xc.xc',
             'run': 'gpaw.cli.main',
             'dos': 'gpaw.gwap.dos',
             'rpa': 'gpaw.xc.rpa.rpa',
             'info': 'gpaw.gwap.info',
             'test': 'gpaw.test.test.main',
             'atom': 'gpaw.atom.aeatom.main',
             'diag': 'gpaw.fulldiag.fulldiag',
             'dataset': 'gpaw.atom.generator2.main',
             'symmetry': 'gpaw.symmetry.analyze_atoms'}


def main():
    commands = functions.keys()
    commands.sort()
    parser1 = optparse.OptionParser(
        usage='Usage: gwap [options] command [more options]\n' +
        '       gwap command --help  (for help on individual commands)',
        description='Run one of these commands: {0}.'
        .format(', '.join(commands)))
    parser1.disable_interspersed_args()
    add = parser1.add_option
    add('-v', '--verbose', action='store_true')
    add('-P', '--parallel', type=int, metavar='N', default=1,
        help="Run on N CPUs.")
    opts1, args1 = parser1.parse_args()

    if opts1.parallel > 1:
        from gpaw.mpi import size
        if size == 1:
            # Start again using gpaw-python in parallel:
            args = ['mpiexec', '-np', str(opts1.parallel),
                    'gpaw-python'] + sys.argv
            os.execvp('mpiexec', args)
    
    if len(args1) == 0:
        parser1.print_help()
        raise SystemExit
    command = args1[0]
    modulename, funcname = functions.get(command, command).rsplit('.', 1)
    module = __import__(modulename, globals(), locals(), [funcname])
    func = getattr(module, funcname)
    kwargs = {}
    if funcname == 'main':
        args = [args1[1:]]
    else:
        nargs, optnames, parser = construct_parser(func, command)
    
        opts, args = parser.parse_args(args1[1:])
        if len(args) != nargs:
            parser.error('Wrong number of arguments!')
        for optname in optnames:
            value = getattr(opts, optname)
            if value is not None:
                kwargs[optname] = value
    try:
        func(*args, **kwargs)
    except Exception as x:
        if opts1.verbose:
            raise
        else:
            print('{0}: {1}'.format(x.__class__.__name__, x), file=sys.stderr)
            print('To get a full traceback, use: gwap --verbose',
                  file=sys.stderr)

            
def construct_parser(func, name):
    """Construct parser from function arguments and docstring."""
    args, varargs, keywords, defaults = inspect.getargspec(func)
    assert varargs is None
    assert keywords is None
    if defaults:
        optnames = args[-len(defaults):]
        defaults = dict(zip(optnames, defaults))
        del args[-len(defaults):]
    else:
        optnames = []
        defaults = {}
    doc = func.func_doc
    headline, doc = doc.split('\n', 1)
    lines = textwrap.dedent(doc).splitlines()
    description = None
    options = []
    shortoptions = set()
    i = 0
    while i < len(lines):
        words = lines[i].split(':')
        if (len(words) > 1 and
            i + 1 < len(lines) and
            lines[i + 1].startswith('    ')):
            if description is None:
                description = ' '.join([headline] + lines[:i])
            arg = words[0]
            help = []
            type = words[1].strip()
            i += 1
            while i < len(lines) and lines[i].startswith('    '):
                help.append(lines[i][4:])
                i += 1
            kwargs = {'help': ' '.join(help)}
            if arg not in defaults:
                continue
            default = defaults.get(arg)
            if type == 'bool':
                kwargs['action'] = 'store_' + str(not default).lower()
            else:
                kwargs['type'] = type
                if default is not None:
                    kwargs['default'] = default
            short = arg[0]
            if short in shortoptions:
                short = arg[1]
            shortoptions.add(short)
            options.append(optparse.Option('-' + short,
                                           '--' + arg.replace('_', '-'),
                                           **kwargs))
        else:
            if description:
                break
            i += 1
        
    epilog = ' '.join(lines[i:])
            
    parser = optparse.OptionParser(
        usage='Usage: gwap {0} <{1}> [options]'.format(name, '> <'.join(args)),
        description=description,
        option_list=options,
        epilog=epilog)
    return len(args), optnames, parser

    
def info(filename):
    """Write summary of GPAW-restart file.
    
    filename: str
        Name of restart-file.
    """
    from gpaw import GPAW
    GPAW(filename)

    
def dos(filename, plot=False, output='dos.csv', width=0.1):
    """Calculate density of states.
    
    filename: str
        Name of restart-file.
    plot: bool
        Show a plot.
    output: str
        Name of CSV output file.
    width: float
        Width of Gaussians.
    """
    from gpaw import GPAW
    from ase.dft.dos import DOS
    calc = GPAW(filename, txt=None)
    dos = DOS(calc, width)
    D = [dos.get_dos(spin) for spin in range(calc.get_number_of_spins())]
    if output:
        fd = sys.stdout if output == '-' else open(output, 'w')
        for x in zip(dos.energies, *D):
            print(*x, sep=', ', file=fd)
        if output != '-':
            fd.close()
    if plot:
        import matplotlib.pyplot as plt
        for y in D:
            plt.plot(dos.energies, y)
        plt.show()
