def agts(queue):
    d = queue.add('diffusion1.py', ncpus=1, walltime=10)
    queue.add('neb.py', deps=d, ncpus=12, walltime=60)
    
