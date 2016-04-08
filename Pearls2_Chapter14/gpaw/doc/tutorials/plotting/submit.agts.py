def agts(queue):
    d = queue.add('CO.py', ncpus=1, walltime=10)
    queue.add('CO2cube.py', deps=d, ncpus=1, walltime=10)
    queue.add('CO2plt.py', deps=d, ncpus=1, walltime=10)
