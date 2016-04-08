def agts(queue):
    calc1 = queue.add('Be_gs_8bands.py',
                      ncpus=2,
                      walltime=20)

    calc2 = queue.add('Be_8bands_lrtddft.py',
                      ncpus=2,
                      walltime=20,
                      deps=calc1)

    calc3 = queue.add('Be_8bands_lrtddft_dE.py',
                      ncpus=2,
                      walltime=20,
                      deps=calc1)

    calc4 = queue.add('Na2_relax_excited.py',
                      ncpus=4,
                      walltime=500)
    
