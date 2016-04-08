import sys
import optparse

from ase.cli.run import Runner, str2dict

from gpaw import GPAW
from gpaw.mixer import Mixer, MixerSum
from gpaw.occupations import FermiDirac, MethfesselPaxton
from gpaw.wavefunctions.pw import PW


class GPAWRunner(Runner):
    def __init__(self):
        Runner.__init__(self)
        self.calculator_name = 'gpaw'

    def make_parser(self):
        parser = optparse.OptionParser(
            usage='gwap run [options] [system, ...]',
            description='Run calculation for system(s).')
        return parser
        
    def add_options(self, parser):
        Runner.add_options(self, parser)
        parser.add_option('-w', '--write', help='Write gpw-file.')
        parser.add_option('-W', '--write-all',
                          help='Write gpw-file with wave functions.')

    def set_calculator(self, atoms, name):
        parameter_namespace = {
            'PW': PW,
            'FermiDirac': FermiDirac,
            'MethfesselPaxton': MethfesselPaxton,
            'Mixer': Mixer,
            'MixerSum': MixerSum}
        parameters = str2dict(self.opts.parameters, parameter_namespace)
        atoms.calc = GPAW(txt=self.get_filename(name, 'txt'), **parameters)

    def calculate(self, atoms, name):
        data = Runner.calculate(self, atoms, name)
        if self.opts.write:
            atoms.calc.write(self.opts.write)
        if self.opts.write_all:
            atoms.calc.write(self.opts.write_all, 'all')
        return data
        
        
def main(args=None):
    runner = GPAWRunner()
    runner.parse(args)
    if runner.errors:
        sys.exit(runner.errors)
