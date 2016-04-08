def agts(queue):
    queue.add('plot_freq.py', creates='nl_freq_grid.png')
    
    simple_si = queue.add('silicon_ABS_simpleversion.py')
    queue.add('plot_silicon_ABS_simple.py', creates='si_abs.png',
              deps=simple_si)
    
    si = queue.add('silicon_ABS.py', creates='mac_eps.csv',
                   ncpus=16, walltime=100)
    queue.add('plot_ABS.py', deps=si, creates='silicon_ABS.png')
    
    al = queue.add('aluminum_EELS.py', ncpus=8, walltime=100)
    queue.add('plot_aluminum_EELS_simple.py', deps=al,
              creates=['aluminum_EELS.png'])
    
    GR = queue.add('graphite_EELS.py', ncpus=8, walltime=100)
    queue.add('plot_EELS.py', deps=GR, creates='graphite_EELS.png')
