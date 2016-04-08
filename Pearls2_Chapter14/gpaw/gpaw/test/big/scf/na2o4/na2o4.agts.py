def agts(queue):
    run = queue.add('na2o4.py', ncpus=4, walltime=2 * 60, deps=[])
