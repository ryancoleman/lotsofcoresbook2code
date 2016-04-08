def agts(queue):
    a = queue.add('LiF_gs.py', ncpus=1, walltime=10)
    queue.add('LiF_RPA.py', ncpus=8, walltime=10, deps=a)
    queue.add('LiF_BSE.py', ncpus=32, walltime=20, deps=a)
