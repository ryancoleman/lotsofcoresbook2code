def agts(queue):
    run = queue.add('C2.py', ncpus=4, walltime=60, deps=[])
