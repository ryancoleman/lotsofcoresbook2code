def agts(queue):
    
    calc1 = queue.add('Na_bulk_gs.py',
                     ncpus=4,
                     walltime=20)

    calc2 = queue.add('Na_bulk_rpa.py',
                     ncpus=16,
                     walltime=500,
                      deps=[calc1])   
 
    queue.add('plot_rpa.py',
              ncpus=1,
              walltime=5,
              deps=[calc1, calc2])
    
 
