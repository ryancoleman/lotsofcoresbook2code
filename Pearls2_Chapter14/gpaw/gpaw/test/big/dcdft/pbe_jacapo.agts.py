def agts(queue):
    run = [queue.add('pbe_jacapo.py %s' % r,
                     queueopts='-l nodes=1:ppn=4:opteron4',
                     ncpus=1,
                     walltime=10*60)
           for r in range(0)]
    if 0:  # do not perform analysis
        # we keep once generated files static
        analyse = queue.add('analyse.py dcdft_pbe_jacapo',
                            ncpus=1, walltime=10, deps=run,
                            creates=['dcdft_pbe_jacapo.csv',
                                     'dcdft_pbe_jacapo.txt',
                                     'dcdft_pbe_jacapo_Delta.txt',
                                     'dcdft_pbe_jacapo_raw.csv'])
