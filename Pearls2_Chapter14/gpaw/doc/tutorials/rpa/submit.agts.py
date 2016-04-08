def agts(queue):
    gs_N2 = queue.add('gs_N2.py', ncpus=1, walltime=2000)
    w = queue.add('frequency.py', deps=gs_N2, walltime=200)
    f = queue.add('con_freq.py', ncpus=2, deps=gs_N2, walltime=1000)
    rpa_N2 = queue.add('rpa_N2.py', deps=gs_N2,
                       #queueopts='-l mem=127GB',  # needed on 16 cpus
                       ncpus=32, walltime=1200)
    queue.add('plot_w.py', deps=w, creates='E_w.png')
    queue.add('plot_con_freq.py', deps=f, creates='con_freq.png')
    queue.add('extrapolate.py', deps=rpa_N2, creates='extrapolate.png')
