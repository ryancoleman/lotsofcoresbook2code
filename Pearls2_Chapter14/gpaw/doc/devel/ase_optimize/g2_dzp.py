import sys

import numpy as np

from ase.data.g2 import data
from ase.data.g2 import molecule_names as names
from ase.tasks.main import run
from ase.tasks.molecule import MoleculeTask

from ase.structure import molecule

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
        m.center()
        return m

class Task(MoleculeTask):
    def __init__(self, **kwargs):
        MoleculeTask.__init__(self, **kwargs)

    def calculate(self, name, atoms):
        data = MoleculeTask.calculate(self, name, atoms)
        data['relaxed fmax'] = np.sqrt((atoms.get_forces()**2).sum(axis=1)).max()
        data['fmax'] = self.fmax
        if 'energy' in data and 'relaxed energy' not in data:
            # no optimization performed
            data['relaxed energy'] = data['energy']
            data['steps'] = 0
        elif 'relaxed energy' in data:
            # optimization
            data['steps'] = self.steps
        if data['relaxed fmax'] > data['fmax']:
            raise RuntimeError('Optimization failed to converge')
        return data

from gpaw.factory import GPAWFactory
from gpaw.mixer import Mixer

if len(sys.argv) == 1:
    optimizer = None
else:
    optimizer = sys.argv[1]

cell = (12.01, 12.02, 12.03)

tag = 'g2_dzp'

if optimizer is not None:
    tag += '_%s' % optimizer

calcopts_default = {
    'mode':'lcao',
    'basis':'dzp',
    'xc':'PBE',
    'width':0.0,
    'fixmom':True,
    'nbands':-2,
    # make systems converge
    'mixer':Mixer(0.05, 2),
    'maxiter':300,
    }
calcfactory_default = GPAWFactory(**calcopts_default)

calcopts_d16 = calcopts_default.copy()
calcopts_d16.update({'convergence': {'density': 1.e-6}})
calcfactory_d16 = GPAWFactory(**calcopts_d16)

calcfactory = calcfactory_default

collection = Collection(data, names, cell)

taskopts = {'fmax': 0.05, 'steps': 100}
if optimizer is not None:
    if 'D16' not in optimizer:
        taskopts.update({'optimizer': optimizer})
        calcfactory = calcfactory_default
    else:
        taskopts.update({'optimizer': optimizer[:-3]})
        calcfactory = calcfactory_d16

task = Task(calcfactory=calcfactory,
            tag=tag,
            use_lock_files=True,
            collection=collection,
            cell=cell,
            **taskopts)
task.run(collection.keys())
