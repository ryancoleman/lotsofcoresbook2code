def agts(queue):
    queue.add('bandstructure.py', ncpus=1, walltime=5,
              creates=['bandstructure.png'])
    
