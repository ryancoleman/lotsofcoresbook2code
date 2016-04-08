def agts(queue):
    # Most of these time out at the moment ...
    return
    
    gs_H_lda = queue.add('H.ralda_01_lda.py', ncpus=2, walltime=5)
    queue.add('H.ralda_02_rpa_at_lda.py', deps=gs_H_lda, ncpus=16, walltime=20)
    queue.add('H.ralda_03_ralda.py', deps=gs_H_lda, ncpus=16, walltime=200)
    gs_H_pbe = queue.add('H.ralda_04_pbe.py', ncpus=2, walltime=5)
    queue.add('H.ralda_05_rpa_at_pbe.py', deps=gs_H_pbe, ncpus=16, walltime=20)
    queue.add('H.ralda_06_rapbe.py', deps=gs_H_pbe, ncpus=16, walltime=200)
    gs_CO = queue.add('CO.ralda_01_pbe+exx.py', ncpus=1, walltime=1000)
    queue.add('CO.ralda_02_CO_rapbe.py', deps=gs_CO, ncpus=16, walltime=2000)
    queue.add('CO.ralda_03_C_rapbe.py', deps=gs_CO, ncpus=16, walltime=2000)
    queue.add('CO.ralda_04_O_rapbe.py', deps=gs_CO, ncpus=16, walltime=2000)
    gs_diamond = queue.add('diamond.ralda_01_pbe.py', deps=gs_CO,
                           ncpus=1, walltime=100)
    queue.add('diamond.ralda_02_rapbe_rpa.py', deps=gs_diamond,
              ncpus=16, walltime=1200)
