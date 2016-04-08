def agts(queue):
    
    calc1 = queue.add('graphite_EELS.py',
                     ncpus=8,
                     walltime=120)

    calc2 = queue.add('silicon_ABS.py',
                     ncpus=16,
                     walltime=40)   
 
    calc3 = queue.add('be_1ml_surf_response.py',
                     ncpus=8,
                     walltime=120)

    queue.add('plot_spectra.py',
              ncpus=1,
              walltime=5,
              deps=[calc1, calc2])
    
    queue.add('check_spectra.py',
              ncpus=1,
              walltime=5,
              deps=[calc1])

