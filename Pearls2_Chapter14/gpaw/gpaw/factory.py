import optparse

from ase.units import Bohr, Hartree
from ase.tasks.calcfactory import CalculatorFactory, str2dict

from gpaw.utilities.gpts import get_number_of_grid_points
from gpaw.wavefunctions.pw import PW


class GPAWFactory(CalculatorFactory):
    def __init__(self, show_text_output=False, write_gpw_file=None,
                 **kwargs):
        self.show_text_output = show_text_output
        # this is used as write(mode=write_gpw_file)
        # so it should be called rather write_gpw_mode?
        if write_gpw_file is not None:
            assert isinstance(write_gpw_file, str)
        self.write_gpw_file = write_gpw_file

        CalculatorFactory.__init__(self, None, 'GPAW', **kwargs)

    def __call__(self, name, atoms):
        kpts = self.calculate_kpts(atoms)

        kwargs = self.kwargs.copy()  # modify a copy
        
        if (not atoms.pbc.any() and len(atoms) == 1 and
            atoms.get_initial_magnetic_moments().any() and
            'hund' not in kwargs):
            kwargs['hund'] = True

        if atoms.pbc.any() and 'gpts' not in kwargs:
            # Use fixed number of gpts:
            h = kwargs.get('h')
            if h is not None:
                h /= Bohr

            cell_cv = atoms.cell / Bohr

            mode = kwargs.get('mode')
            if mode == 'pw':
                mode = PW()

            gpts = get_number_of_grid_points(cell_cv, h, mode,
                                             kwargs.get('realspace'))
            kwargs['h'] = None
            kwargs['gpts'] = gpts
            
            if isinstance(mode, PW):
                kwargs['mode'] = PW(mode.ecut * Hartree,
                                    mode.fftwflags,
                                    atoms.cell)

        if self.show_text_output:
            txt = '-'
        else:
            txt = name + '.txt'


        from gpaw import GPAW
        return GPAW(txt=txt, kpts=kpts, **kwargs)
        
    def add_options(self, parser):
        CalculatorFactory.add_options(self, parser)
        
        calc = optparse.OptionGroup(parser, 'GPAW')
        calc.add_option('--parameter-file', metavar='FILE',
                        help='Read GPAW parameters from file.')
        calc.add_option('-S', '--show-text-output', action='store_true',
                        help='Send text output from calculation to ' +
                        'standard out.')
        calc.add_option('-W', '--write-gpw-file', metavar='MODE',
                        help='Write gpw file.')
        parser.add_option_group(calc)

    def parse(self, opts, args):
        if opts.parameters:
            # Import stuff that eval() may need to know:
            from gpaw.wavefunctions.pw import PW
            from gpaw.occupations import FermiDirac, MethfesselPaxton
            from gpaw.mixer import Mixer, MixerSum
            from gpaw.poisson import PoissonSolver
            from gpaw.eigensolvers import RMM_DIIS
       
            self.kwargs.update(str2dict(opts.parameters, locals()))
            opts.parameters = None

        CalculatorFactory.parse(self, opts, args)

        self.show_text_output = opts.show_text_output
        self.write_gpw_file = opts.write_gpw_file

        if opts.parameter_file:
            self.kwargs.update(eval(open(opts.parameter_file).read()))
