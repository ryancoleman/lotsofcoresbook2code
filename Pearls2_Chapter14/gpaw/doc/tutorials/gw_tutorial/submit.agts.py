def agts(queue):
    groundstate = queue.add('Si_groundstate.py')
    queue.add('Si_gw.py', deps=groundstate, ncpus=1, walltime=20)
    conv = queue.add('exx_convergence.py', ncpus=4, walltime=60)
    queue.add('plot_convergence.py', deps=conv, creates=['Si_EXX.png'])
    freq = queue.add('frequency.py', ncpus=4, walltime=40*60)
    queue.add('plot_frequency.py', deps=freq, creates='Si_w.png')
