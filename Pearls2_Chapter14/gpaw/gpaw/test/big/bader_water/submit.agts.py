def agts(queue):
    calc = queue.add('bader_water.py', ncpus=8, walltime=10)
    queue.add('bader_plot.py', walltime=5, deps=[calc])

