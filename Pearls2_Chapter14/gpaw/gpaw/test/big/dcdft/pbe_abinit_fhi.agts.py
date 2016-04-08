def agts(queue):
    run = [queue.add('pbe_abinit_fhi.py %s' % r,
                     queueopts='-l nodes=1:ppn=4:opteron4',
                     ncpus=1,
                     walltime=5*60)
           for r in range(1)]
    if 0:  # do not perform analysis
        # we keep once generated files static
        analyse = queue.add('analyse.py dcdft_pbe_abinit_fhi',
                            ncpus=1, walltime=10, deps=run,
                            creates=['dcdft_pbe_abinit_fhi.csv',
                                     'dcdft_pbe_abinit_fhi.txt',
                                     'dcdft_pbe_abinit_fhi_Delta.txt',
                                     'dcdft_pbe_abinit_fhi_raw.csv'])
