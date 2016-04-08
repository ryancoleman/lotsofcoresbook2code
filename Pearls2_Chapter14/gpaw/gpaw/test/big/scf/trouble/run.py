import os
import optparse
import traceback
from time import time
from glob import glob

import ase.db

import gpaw.mpi
from gpaw import GPAW


op = optparse.OptionParser(
    usage='python run.py dir [tes1.py test2.py ...]',
    description='Run tests in dir which must contain a file called ' +
    'params.py defining a calc(atoms) function that sets the eigensolver ' +
    'and other stuff')
opts, args = op.parse_args()

tag = args.pop(0).rstrip('/')
execfile(tag + '/params.py')  # get the calc function
        
c = ase.db.connect('results.db')

if args:
    names = args
else:
    names = [name for name in glob('*.py') if name not in
             ['run.py', 'params.py',
              'submit.agts.py', 'run.py.py', 'analyse.py']]
    
for name in names:
    namespace = {}
    execfile(name, namespace)
    atoms = namespace['atoms']
    ncpus = namespace.get('ncpus', 8)
    
    if gpaw.mpi.size != ncpus:
        continue
        
    name = name[:-3]

    id = c.reserve(name=name, tag=tag)
    if not id:
        continue
        
    if atoms.calc is None:
        atoms.calc = GPAW()

    calc(atoms)
    
    atoms.calc.set(txt=tag + '/' + name + '.txt')
    
    t = time()
    try:
        e1 = atoms.get_potential_energy()
        ok = True
    except:
        ok = False
        if gpaw.mpi.rank == 0:
            traceback.print_exc(file=open(tag + '/' + name + '.error', 'w'))
    t = time() - t

    c.write(atoms, name=name, tag=tag, ok=ok,
            time=t, iters=atoms.calc.iter, ncpus=ncpus)
    
    del c[id]
