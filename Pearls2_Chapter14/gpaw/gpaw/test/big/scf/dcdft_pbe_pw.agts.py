runs = [
#    'bands00m',  # the default number of bands is bad
#    #'m051',  # slower mixing does not help
#    #'m103poisson',  # no influence
#    'm103', 'm105', 'm107',  # for Mixer nmaxold matters: large is bad
#    'b103', 'd103', 's103',
#    #'d105', 'd107', 'd203',  # for MixerDif nmaxold does not matter
#    'inititer10m', # cg steps and switching to rmm-diis does not help
#    'cgm103', 'cgm203',  # cg helps
#    #'cgs103',  # cg + MixerSum
#    'cgd103', 'cgd203', 'cgd253', # cg + MixerDif: the winner
#    'dzpm103',  # lcao guess does not help
#    #'cgdzpm103',  # lcao + cg
#    #'cgdzps103',  # lcao + cg + MixerSum
#    #'cgdzpd103',  # lcao + cg + MixerDif
#    'fm103',  # FFT
    'davd203',  # PW mode winner
    'davfd203',  # FFT
    ]
runsstr = ','.join(runs)
runsstr += ',jacapo'

def agts(queue):
    run = [queue.add('dcdft_pbe_pw.py %s --gpaw=fprojectors=1' % r,
                     ncpus=1,
                     walltime=30*60)
           for r in runs*8]
    jacapo = queue.add('dcdft_pbe_pw_jacapo.py',
                       ncpus=1,
                       walltime=10*60)
    if 0:  # we don't run all the runs for the moment
        analyse = queue.add('analyse.py bulk scf_dcdft_pbe_pw ' + runsstr,
                            ncpus=1, walltime=10, deps=run + [jacapo],
                            creates=['scf_dcdft_pbe_pw_energy.csv',
                                     'scf_dcdft_pbe_pw_calculator_steps.png'])
