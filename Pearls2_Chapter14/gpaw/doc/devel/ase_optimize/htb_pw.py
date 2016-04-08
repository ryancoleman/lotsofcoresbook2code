import sys

import numpy as np

from ase.test.tasks.htb import HaasTranBlahaBulkTask

class Task(HaasTranBlahaBulkTask):
    def __init__(self, xc='LDA', **kwargs):
        HaasTranBlahaBulkTask.__init__(self, xc, **kwargs)

    def calculate(self, name, atoms):
        data = HaasTranBlahaBulkTask.calculate(self, name, atoms)
        sforces = - atoms.get_stress().ravel()*atoms.get_volume()
        rmssforces = np.sum(sforces**2)**0.5
        data['relaxed fmax'] = rmssforces
        data['fmax'] = self.sfmax
        if 'energy' in data and 'relaxed energy' not in data:
            # no optimization performed
            data['relaxed energy'] = data['energy']
            data['steps'] = 0
        elif 'relaxed energy' in data:
            # optimization
            data['steps'] = self.ssteps
        if 'strain optimizer force calls' in data:
            data['optimizer force calls'] = data['strain optimizer force calls']
        if data['relaxed fmax'] > data['fmax']:
            raise RuntimeError('Optimization failed to converge')
        return data

from gpaw import PW
from gpaw.mixer import Mixer
from gpaw.factory import GPAWFactory
from gpaw.mpi import serial_comm

if len(sys.argv) == 1:
    optimizer = None
else:
    optimizer = sys.argv[1]

tag = 'htb_pw'

if optimizer is not None:
    tag += '_%s' % optimizer

calcfactory = GPAWFactory(
    mode=PW(),
    xc='PBE',
    width=0.1,
    maxiter=400,
    mixer=Mixer(0.10, 2),
    eigensolver='cg',
    # avoid problems with band parallelization
    communicator=serial_comm,
    )

taskopts = {'sfmax': 0.01, 'ssteps': 50}
if optimizer is not None:
    taskopts.update({'soptimizer': optimizer})

task = Task(xc='LDA',
            calcfactory=calcfactory,
            tag=tag,
            use_lock_files=True,
            **taskopts)

keys = task.collection.names
# we don't have setups for those
for s in ['CeO2', 'Th']:
    keys.remove(s)
task.run(keys)
