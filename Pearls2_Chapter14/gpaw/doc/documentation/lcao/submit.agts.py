def agts(queue):
    queue.add('basisgeneration.py', ncpus=1, walltime=10)
    queue.add('lcao_h2o.py', ncpus=1, walltime=10)
    queue.add('lcao_opt.py', ncpus=1, walltime=10)

