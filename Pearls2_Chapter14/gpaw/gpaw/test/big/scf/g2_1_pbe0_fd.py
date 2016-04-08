import sys

from ase.data.g2_1 import data
from ase.data.g2_1 import atom_names
from ase.data.g2_1 import molecule_names
from ase.tasks.main import run
from ase.tasks.molecule import MoleculeTask

from ase.structure import molecule

from gpaw.test.big.scf.analyse import rundefs

class Collection:

    def __init__(self, data, names, cell):
        self.data = data
        self.names = names
        self.cell = cell

    def __getitem__(self, name):
        return self.create_item(name)

    def keys(self):
        return self.names

    def create_item(self, name):
        m = molecule(name, data=self.data)
        m.set_cell(self.cell)
        m.set_pbc((0, 0, 0))
        m.center()
        return m

class Task(MoleculeTask):
    def __init__(self, **kwargs):
        MoleculeTask.__init__(self, **kwargs)

    def calculate(self, name, atoms):
        data = MoleculeTask.calculate(self, name, atoms)
        try:
            steps = atoms.get_calculator().get_number_of_iterations()
        except (AttributeError, NotImplemented):
            steps = None
        data['calculator steps'] = steps
        return data

from gpaw.mixer import Mixer, MixerSum, MixerSum2, MixerDif
from gpaw.mixer import BroydenMixer, BroydenMixerSum
from gpaw.factory import GPAWFactory

if len(sys.argv) == 1:
    run = None
else:
    run = sys.argv[1]

cell = (12.01, 12.02, 12.03)

tag = 'scf_g2_1_pbe0_fd'

if run is not None:
    tag += '_%s' % run

calcopts = {
    'mode': 'fd',
    'xc': 'PBE0',
    'width': 0.0,
    'fixmom': True,
    # allow other mixers
    'spinpol': True,
    # allow for long SCFs
    'maxiter': 700,
    }

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
else:
    calcopts.update({'eigensolver': 'rmm-diis'})
if run.startswith('cgdzp'):
    calcopts.update({'basis': 'dzp'})
calcopts.update({'mixer': eval(rundefs[run])})

calcfactory = GPAWFactory(**calcopts)

collection = Collection(data, atom_names + molecule_names, cell)

task = Task(calcfactory=calcfactory,
            tag=tag,
            use_lock_files=True,
            collection=collection,
            cell=cell,
            )
keys = collection.keys()
for m in ['Na2', 'NaCl']:  # those seem to need cg
    keys.remove(m)
task.run(keys)
