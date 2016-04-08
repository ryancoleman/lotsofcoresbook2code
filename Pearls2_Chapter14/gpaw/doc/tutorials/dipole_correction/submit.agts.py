def agts(queue):
    d = queue.add('dipole.py', ncpus=4, walltime=60)
    queue.add('plot.py', deps=d, ncpus=1, walltime=10,
              creates=['zero.png', 'periodic.png', 'corrected.png',
                       'slab.png'])
    queue.add('check.py', deps=d, ncpus=1, walltime=10)
