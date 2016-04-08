def agts(queue):
   
    calc = queue.add('Si_AH.py', walltime=120, ncpus=1)
    queue.add('Si_bandstructure.py', walltime=12 * 60, ncpus=8, deps=calc)
