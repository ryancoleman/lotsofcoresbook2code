def agts(queue):
    queue.add('pbe_aims.py Al', ncpus=1,
              queueopts='-l nodes=1:ppn=1:xeon16', walltime=40)
