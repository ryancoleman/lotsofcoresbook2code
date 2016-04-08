def agts(queue):
    bulk = queue.add('bulk.py', ncpus=4, walltime=6)
    surf = queue.add('surface.py', ncpus=4, walltime=6)
    sigma = queue.add('sigma.py', deps=[bulk, surf])
    queue.add('fig2.py', deps=sigma, creates='fig2.png')
