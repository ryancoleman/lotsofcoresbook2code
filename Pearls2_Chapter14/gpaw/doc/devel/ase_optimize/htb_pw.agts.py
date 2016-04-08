optimizers = [
    'BFGS',
    'LBFGS',
    'FIRE',
    #'MDMin',  # memory hungy!
    #'BFGSLineSearch',  # StrainFilter instance has no attribute 'copy'
    #'LBFGSLineSearch',  # hundreds of force evaluations
    #'SciPyFminCG',  # hundreds of force evaluations
    #'SciPyFminBFGS',  # scipy 0.12.0: hundreds of force evaluations
    'GoodOldQuasiNewton',
]

runsstr = ','.join(optimizers)

def agts(queue):
    run = [queue.add('htb_pw.py %s' % o,
                     ncpus=1,
                     walltime=40*60)
           for o in optimizers]
    analyse = queue.add('task_analyse.py bulk htb_pw ' + runsstr,
                        ncpus=1, walltime=10, deps=run,
                        creates=['htb_pw_relaxed_energy.csv',
                                 'htb_pw_optimizer_force_calls.png'])
