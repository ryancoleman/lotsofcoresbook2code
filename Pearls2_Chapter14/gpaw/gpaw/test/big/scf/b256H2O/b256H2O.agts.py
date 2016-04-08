def agts(queue):
    runA = queue.add('b256H2O.py A', ncpus=4, walltime=5*60, deps=[])
    runB = queue.add('b256H2O.py B', ncpus=4, walltime=5*60, deps=[])
