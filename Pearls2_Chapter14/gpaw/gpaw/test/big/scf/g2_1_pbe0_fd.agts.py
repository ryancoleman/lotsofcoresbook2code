runs = [
    #'bands05m', 'bands10m',  # 'bands20m',  # large number of bands bad
    #'m101',
    #'m102',
    #'m103', 'm105',  # 'm107',  # for Mixer nmaxold matters
    #'m051',  # many steps, and no advantage for this set of systems
    #'m203',  # larger mixing better
    's103',
    'd103', 'd203', 'd253',  # MixerDiff best
    'dzpm103',  'dzpm203',  # 'dzpm253', # dzp guess does not help
    ]
runsstr = ','.join(runs)

def agts(queue):
    run = [queue.add('g2_1_pbe0_fd.py %s --gpaw=fprojectors=1' % r,
                     ncpus=4,
                     walltime=40*60)
           for r in runs * 2]
    analyse = queue.add('analyse.py molecule scf_g2_1_pbe0_fd ' + runsstr,
                        ncpus=1, walltime=10, deps=run,
                        creates=['scf_g2_1_pbe0_fd_energy.csv',
                                 'scf_g2_1_pbe0_fd_calculator_steps.png'])
