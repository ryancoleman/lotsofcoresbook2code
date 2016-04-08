import sys

from ase.test.tasks.dcdft import DeltaCodesDFTTask

from gpaw.atom.configurations import parameters
from gpaw.test.big.scf.analyse import rundefs

if len(sys.argv) == 1:
    xc = 'PBE'
    run = None
elif len(sys.argv) == 2:
    xc = 'PBE'
    run = sys.argv[1]
elif len(sys.argv) == 3:
    xc = sys.argv[1]
    run = sys.argv[2]

tag = 'scf_dcdft_%s_pw' % (xc.lower())

if run is not None:
    tag += '_%s' % run

class Task(DeltaCodesDFTTask):
    def __init__(self, inititer=0, **kwargs):
        DeltaCodesDFTTask.__init__(self, **kwargs)
        self.inititer = inititer

    def calculate(self, name, atoms):
        data = DeltaCodesDFTTask.calculate(self, name, atoms)
        try:
            steps = atoms.get_calculator().get_number_of_iterations()
        except (AttributeError, NotImplemented):
            steps = None
        data['calculator steps'] = steps + self.inititer
        return data

from gpaw import ConvergenceError
from gpaw import PW
from gpaw.factory import GPAWFactory
from gpaw.mixer import Mixer, MixerSum, MixerSum2, MixerDif
from gpaw.mixer import FFTMixer, FFTMixerSum, FFTMixerDif
from gpaw.mixer import BroydenMixer, BroydenMixerSum
from gpaw.poisson import PoissonSolver
from gpaw.occupations import FermiDirac, MethfesselPaxton

class Factory(GPAWFactory):
    def __init__(self, show_text_output=False, write_gpw_file=None,
                 inititer=0, **kwargs):
        GPAWFactory.__init__(self, show_text_output=show_text_output,
                             write_gpw_file=write_gpw_file,
                             **kwargs)
        self.inititer = inititer
        self.maxiter = kwargs['maxiter']
        self.eigensolver = kwargs['eigensolver']

    def __call__(self, name, atoms):
        calculator = GPAWFactory.__call__(self, name, atoms)
        if name.split('-')[0] in ['Li', 'Na']:
            # https://listserv.fysik.dtu.dk/pipermail/gpaw-developers/2012-May/002870.html
            calculator.set(h=0.11)
        if self.inititer > 0:
            try:
                calculator.set(eigensolver='cg')
                calculator.set(maxiter=self.inititer)
                atoms.set_calculator(calculator)
                atoms.get_potential_energy()
            except ConvergenceError:
                pass
            calculator.set(maxiter=self.maxiter)
            calculator.set(eigensolver=self.eigensolver)
        return calculator

calcopts = {
    'mode': PW(),
    'xc': xc,
    # allow other mixers
    'spinpol': True,
    # allow for long SCFs
    'maxiter': 500,
    'nbands': -5,
    }

if run.startswith('inititer'):
    inititer = int(run[len('inititer'):len('inititer') + 2])
    calcopts.update({'inititer': inititer})
else:
    inititer = 0
if run.startswith('bands'):
    nbands = run[len('bands'):len('bands') + 2]
    if nbands == '00':
        calcopts.update({'nbands': None})
    else:
        calcopts.update({'nbands': - int(nbands)})
if run.startswith('cgbands'):
    nbands = run[len('cgbands'):len('cgbands') + 2]
    if nbands == '00':
        calcopts.update({'nbands': None})
    else:
        calcopts.update({'nbands': - int(nbands)})
if run.startswith('dzpbands'):
    nbands = run[len('dzpbands'):len('dzpbands') + 2]
    if nbands == '00':
        calcopts.update({'nbands': None})
    else:
        calcopts.update({'nbands': - int(nbands)})
if run.startswith('szpdzp'):
    calcopts.update({'basis': 'szp(dzp)'})
if run.startswith('szdzp'):
    calcopts.update({'basis': 'sz(dzp)'})
if run.startswith('dzp'):
    calcopts.update({'basis': 'dzp'})
if 'cg' in run:
    calcopts.update({'eigensolver': 'cg'})
elif 'dav' in run:
    calcopts.update({'eigensolver': 'dav'})
else:
    calcopts.update({'eigensolver': 'rmm-diis'})
if run.startswith('cgdzp'):
    calcopts.update({'basis': 'dzp'})
calcopts.update({'mixer': eval(rundefs[run])})
if 'mp' in run:
    calcopts.update({'occupations': MethfesselPaxton(0.1)})
if 'poisson' in run:
    calcopts.update({'poissonsolver': PoissonSolver(eps=1e-12)})

calcfactory = Factory(**calcopts)

task = Task(
    inititer=inititer,
    calcfactory=calcfactory,
    tag=tag,
    use_lock_files=True,
    )

if __name__ == '__main__':
    # run systems from collection for which we have setups
    keys = list(set(parameters.keys()).intersection(task.collection.names))
    keys.sort()
    task.run(keys)
